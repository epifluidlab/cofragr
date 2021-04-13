#' @export
read_fragments <- function(file_path, range = NULL, genome = NULL) {
  frag <- bedtorch::read_bed(file_path = file_path,
                             range = range,
                             genome = genome)
  frag_metadata <- mcols(frag)

  mapq_guessed <- FALSE
  if (!"mapq" %in% colnames(frag_metadata)) {
    logging::logwarn("mapq not in column names. Need to guess which column is mapq")

    for (col_idx in seq_along(colnames(frag_metadata))) {
      if (is.integer(frag_metadata[[col_idx]])) {
        logging::logwarn(str_interp("Column ${colnames(frag_metadata)[col_idx]} seems to be mapq"))
        colnames(frag_metadata)[col_idx] <- "mapq"
        mapq_guessed <- TRUE
        break
      }
    }

    if (!"mapq" %in% colnames(frag_metadata))
      stop("Cannot find mapq in the data frame")
  }

  # Only need mapq. Drop other columns
  mcols(frag) <- frag_metadata["mapq"]
  if (mapq_guessed)
    print(frag)
  frag
}

#' @export
write_contact_matrix <- function(cm, file_path, comments = NULL) {
  cm <- as_tibble(cm) %>%
    select(-c(width, strand, ranges2.width)) %>%
    rename(chrom = seqnames, chrom2 = seqnames2, start2 = ranges2.start, end2 = ranges2.end) %>%
    select(c(chrom, start, end, chrom2, start2, end2, everything())) %>%
    mutate(start = as.integer(start - 1), start2 = as.integer(start2 - 1))

  cm %>% data.table::as.data.table() %>% bedtorch::write_bed(file_path, comments = comments)
}

# write_contact_matrix(tmp3) %>% bedtorch::write_bed(file_path = "cm.bed")

# Aggregate fragment lengths over the window
# Each row is a genomic window, and len is a vector of fragment lengths
# Fragment-over-the-window is determined by its mid-point
# # A tibble: 12 x 2
# start len
# <int> <list>
#     1 29900000 <int [1]>
#     2 30000000 <int [1,265]>
#     3 30100000 <int [1,313]>
fraglen_over_window <- function(frag, bin_size = 500e3L) {
  stopifnot(bin_size %in% c(100e3L, 500e3L, 1e6L, 2.5e6L, 5e6L))

  # Only allow single chromosome fragment data
  chrom <- unique(seqnames(frag))
  stopifnot(length(chrom) == 1)

  # if (is.null(genome)) {
  #   # Get the genome name. E.g. GRCh37
  #   genome <- unique(GenomeInfoDb::genome(frag))
  #   stopifnot(length(genome) == 1)
  #   if (!is.na(genome))
  #     genome <- as.character(genome)
  # }

  # Find mid-points
  midpoints <-
    as.integer(round(start(frag) + (width(frag) - 1) / 2))
  frag$len <- width(frag)
  start(frag) <- midpoints
  width(frag) <- 1

  # if (is.null(genome) || is.na(genome)) {
    window <-
      bedtorch::make_windows(
        window_size = bin_size,
        chrom = chrom,
        chrom_sizes = data.table::data.table(chrom = as.character(chrom),
                                             size = as.integer(max(end(
                                               frag
                                             ))))
      )
  # }
  # else {
  #   window <-
  #     bedtorch::make_windows(window_size = bin_size,
  #                            chrom = chrom,
  #                            genome = genome)
  # }

  hits <- GenomicRanges::findOverlaps(frag, window)
  as_tibble(hits) %>%
    mutate(len = frag$len[queryHits]) %>%
    group_by(subjectHits) %>%
    summarize(len = list(len)) %>%
    mutate(start = as.integer(start(window)[subjectHits] - 1)) %>%
    select(start, len)
}


# Build the calculation plan: each row is a pair of bins, with the indices to
# look up for fraglen
# # A tibble: 2,556 x 4
#      start1   start2  id.1  id.2
#      <int>    <int> <int> <int>
#  1 16000000 16000000     1     1
#  2 16000000 16500000     1     2
#  3 16000000 17000000     1     3
#  4 16000000 17500000     1     4
cofrag_plan <- function(fraglen, bin_size = 500e3L) {
  min_start <- min(fraglen$start)
  max_start <- max(fraglen$start)

  stopifnot(max_start > min_start && (max_start - min_start) %% bin_size == 0)

  plan <- expand_grid(start1 = seq(min_start, max_start, by = bin_size),
                      start2 = seq(min_start, max_start, by = bin_size)) %>%
    filter(start1 < start2) %>%
    arrange(start1, start2)

  fraglen %<>% select(start) %>% mutate(id = seq.int(start))

  plan %>%
    inner_join(fraglen, by = c(start1 = "start")) %>%
    rename(id.1 = id) %>%
    inner_join(fraglen, by = c(start2 = "start")) %>%
    rename(id.2 = id)
}


stat_func_ks <- function(subsample, min_sample_size = NULL, bootstrap = 1L) {
  function(len1, len2) {
    stopifnot(bootstrap >= 1)

    n_frag1 <- length(len1)
    n_frag2 <- length(len2)
    if (min(n_frag1, n_frag2) < min_sample_size)
      return(NULL)
      # return(tibble(n_frag1 = n_frag1, n_frag2 = n_frag2, bootstrap = NA, p_value = NA))

    v <- seq.int(bootstrap) %>% map_dbl(function(idx) {
      len1 <- sample(len1, size = subsample, replace = TRUE)
      len2 <- sample(len2, size = subsample, replace = TRUE)
      suppressWarnings(ks.test(len1, len2)$p.value)
    })
    p_value <- mean(v)
    p_value_sd <- sd(v)

    tibble(
      score = min(16,-log10(p_value)),
      n_frag1 = n_frag1,
      n_frag2 = n_frag2,
      p_value = p_value,
      p_value_sd = p_value_sd
    )
  }
}


calc_contact_matrix <- function(fraglen, cofrag_plan, stat_func, bin_size = 500e3L, n_workers = 1L, seed = NULL) {
  if (n_workers > 1) {
    cl <- makeForkCluster(n_workers)
    registerDoParallel(cl)
    logging::loginfo(str_interp("Successfully registered a cluster with ${n_workers} workers"))
    on.exit(stopCluster(cl), add = TRUE)
  }

  # Divide the plan among workers
  plan_idx = unique(floor(seq.int(1, nrow(cofrag_plan) + 1, length.out = n_workers + 1)))
  # Sanity check
  stopifnot(all(tail(plan_idx, -1) > tail(lag(plan_idx), -1)))

  if (!is.null(seed))
    set.seed(seed)

  foreach(
    idx_start = head(plan_idx,-1),
    idx_end = tail(plan_idx,-1) - 1,
    .combine = "rbind",
    .multicombine = TRUE
    # .export = c("fraglen", "cofrag_plan")
  ) %dorng% {
    id_list_1 <- cofrag_plan$id.1[idx_start:idx_end]
    id_list_2 <- cofrag_plan$id.2[idx_start:idx_end]

    map2_dfr(id_list_1,
             id_list_2,
             function(id.1, id.2) {
               len1 <- fraglen$len[[id.1]]
               len2 <- fraglen$len[[id.2]]
               result <- stat_func(len1, len2)
               if (!is.null(result))
                 result %<>% mutate(id.1 = id.1, id.2 = id.2)
               result
             })
  } %>%
    inner_join(x = cofrag_plan, y = ., by = c("id.1", "id.2")) %>%
    select(-c(id.1, id.2))
    ## Cap the score at 16
    # mutate(score0 = pmin(16, -log10(p_value))) %>%
    ## Aggregate over bootstrap iterations
    # group_by(start1, start2) %>%
    # mutate(score = pmin(16, -log10(mean(p_value)))) %>%
    # ungroup()
}


build_gr <- function(contact_matrix, chrom, bin_size, seqinfo = NULL) {
  stopifnot(length(chrom) == 1)

  gr <- GenomicRanges::GRanges(seqnames = chrom,
                               IRanges::IRanges(contact_matrix$start1 + 1, width = bin_size),
                               seqinfo = seqinfo)

  ranges2 <- IRanges::IRanges(contact_matrix$start2 + 1, width = bin_size)
  df <- contact_matrix%>% select(-c(1:2)) %>% DataFrame()
  df$n_frag1 <- Rle(df$n_frag1)
  df$n_frag2 <- Rle(df$n_frag2)
  mcols(gr) <- cbind(DataFrame(seqnames2 = Rle(chrom, nrow(contact_matrix))), ranges2, df)
  gr
}

#' @export
contact_matrix <-
  function(frag,
           bin_size = 500e3L,
           n_workers = 1L,
           subsample = 10e3L,
           min_sample_size = 100L,
           bootstrap = 1L,
           seed = NULL) {
    chrom <- unique(seqnames(frag))
    stopifnot(length(chrom) >= 1)

    chrom %>% map(function(chrom) {
      frag <- frag[seqnames(frag) == chrom]

      fraglen <-
        fraglen_over_window(frag, bin_size = bin_size)

      plan <- cofrag_plan(fraglen, bin_size = bin_size)
      calc_contact_matrix(
        fraglen,
        plan,
        n_workers = n_workers,
        stat_func = stat_func_ks(
          subsample = subsample,
          min_sample_size = min_sample_size,
          bootstrap = bootstrap
        ),
        bin_size = bin_size,
        seed = seed
      ) %>%
        build_gr(chrom = chrom, bin_size = bin_size, seqinfo = GenomicRanges::seqinfo(frag))
    }) %>%
      do.call(c, args = .)
  }
