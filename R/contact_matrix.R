#' Read a fragment data file
#'
#' A fragment data file is essentially just a BED file. The reason we need a
#' special function here to load a fragment data file is that, by convention,
#' the data frame should contain a column representing MAPQ scores. However, the
#' data file may not have a column header. Therefore, we may need to guess which
#' column represents MAPQ scores.
#'
#' @param file_path path to the file.
#' @param range read the entire file if `range` is `NULL`.
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

  # Find mid-points
  midpoints <-
    as.integer(round(start(frag) + (width(frag) - 1) / 2))
  frag$len <- width(frag)
  start(frag) <- midpoints
  width(frag) <- 1

  window <-
    bedtorch::make_windows(
      window_size = bin_size,
      chrom = chrom,
      genome = GenomeInfoDb::genome(seqinfo(frag)) %>% unique()
    )

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
    filter(start1 <= start2) %>%
    arrange(start1, start2)

  fraglen %<>% select(start) %>% mutate(id = seq.int(start))

  plan %>%
    inner_join(fraglen, by = c(start1 = "start")) %>%
    rename(id.1 = id) %>%
    inner_join(fraglen, by = c(start2 = "start")) %>%
    rename(id.2 = id)
}


# A function factory that generates a stat_func, which is used to calculate the
# distance metrics between two collections of fragments
stat_func_ks <- function(subsample, min_sample_size = NULL, bootstrap = 1L) {
  # We need to force early-evaluation, since the stat_func may be used in a
  # %dopar% worker context. Failed to do so may leads to "object can't be found"
  # error.
  force(subsample)
  force(min_sample_size)
  force(bootstrap)

  function(len1, len2) {
    stopifnot(bootstrap >= 1)

    n_frag1 <- length(len1)
    n_frag2 <- length(len2)
    if (min(n_frag1, n_frag2) < min_sample_size)
      return(NULL)
      # return(tibble(n_frag1 = n_frag1, n_frag2 = n_frag2, bootstrap = NA, p_value = NA))

    v <- purrr::map_dbl(seq.int(bootstrap), function(idx) {
      # v <- seq.int(bootstrap) %>% map_dbl(function(idx) {
      len1 <- sample(len1, size = subsample, replace = TRUE)
      len2 <- sample(len2, size = subsample, replace = TRUE)
      suppressWarnings(ks.test(len1, len2)$p.value)
    })
    # Attention: if you calculate mean(pvalues) and then calculate score, it doesn't work
    # The compartment scores will look weird.
    score <- median(-log10(2e-16) + log10(pmax(2e-16, v)))
    p_value <- mean(v)
    p_value_sd <- sd(v)

    dplyr::tibble(
      score = score,
      n_frag1 = n_frag1,
      n_frag2 = n_frag2,
      p_value = p_value,
      p_value_sd = p_value_sd
    )
  }
}


calc_contact_matrix <- function(fraglen,
                                cofrag_plan,
                                stat_func,
                                bin_size = 500e3L,
                                block_size = 10e6L,
                                n_workers = 1L,
                                seed = NULL) {
  if (n_workers > 1) {
    # cl <- makeForkCluster(n_workers, outfile = "")
    cl <- makeCluster(n_workers, outfile = "")
    registerDoParallel(cl)
    logging::loginfo(str_interp("Successfully registered a cluster with ${n_workers} workers"))
    on.exit(stopCluster(cl), add = TRUE)
  }

  # Divide fraglen into blocks
  n_blocks <- (max(fraglen$start) - min(fraglen$start) + bin_size) / block_size
  n_blocks <- max(as.integer(round(n_blocks * (n_blocks + 1) / 2)), n_workers)
  logging::loginfo(str_interp("Divided the task into ${n_blocks} segments"))

  # Divide the cofrag_plan among workers
  # plan_idx are the start and end boundaries of the segments
  plan_idx = unique(floor(seq.int(1, nrow(cofrag_plan) + 1, length.out = n_blocks + 1)))
  # Sanity check
  stopifnot(all(tail(plan_idx,-1) > tail(lag(plan_idx),-1)))

  # This is the main data input for doParallel jobs
  # Each object in foreach_layout is a segment of cofrag_plan, and is provided
  # to a particular worker
  foreach_layout <- 1:(length(plan_idx) - 1) %>%
    map(function(idx) {
      id_list_1 <- cofrag_plan$id.1[plan_idx[idx]:(plan_idx[idx + 1] - 1)]
      id_list_2 <-
        cofrag_plan$id.2[plan_idx[idx]:(plan_idx[idx + 1] - 1)]
      # Only keep fraglen rows that will be used in the job
      fraglen_map <- new.env()
      unique(c(id_list_1, id_list_2)) %>%
        walk(function(id) {
          fraglen_map[[as.character(id)]] <- fraglen$len[[id]]
        })

      list(
        block_id = idx,
        id_list_1 = id_list_1,
        id_list_2 = id_list_2,
        fraglen_map = fraglen_map
      )
    })

  if (!is.null(seed))
    set.seed(seed)

  foreach(
    data = foreach_layout,
    .export = "n_blocks",
    .combine = "rbind",
    .multicombine = TRUE
  ) %dorng% {
    block_id <- data$block_id

    worker_name <-
      paste(Sys.info()[['nodename']], Sys.getpid(), sep = '-')
    mem_mb <- as.numeric(lobstr::mem_used()) / 1e6
    logging::loginfo(stringr::str_interp("${worker_name}: processing ${block_id}/${n_blocks}, memory used: ${mem_mb} MB"))

    id_list_1 <- data$id_list_1
    id_list_2 <- data$id_list_2
    fraglen_map <- data$fraglen_map

    stopifnot(length(id_list_1) == length(id_list_2))

    purrr::map_dfr(1:length(id_list_1), function(idx) {
      id.1 <- id_list_1[idx]
      id.2 <- id_list_2[idx]

      len1 <- fraglen_map[[as.character(id.1)]]
      len2 <- fraglen_map[[as.character(id.2)]]

      result <- stat_func(len1, len2)
      if (!is.null(result))
        result <-
        dplyr::mutate(result, id.1 = id.1, id.2 = id.2)

      result
    })
  } %>%
    inner_join(x = cofrag_plan,
               y = .,
               by = c("id.1", "id.2")) %>%
    select(-c(id.1, id.2))
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


#' Preprocess the fragment BED-like data frame
#'
#' This function takes a fragment BED as input. For each chromosome, return a
#' new data frame where each row is for a specific position, associated with a
#' vector that represents the lengths of all fragments at that position. The
#' result is much more compact than the original BED-like data frame.
#'
#' Example of returned value:

#' $`21`
#' # A tibble: 74 x 2
#' start len
#' <int> <list>
#'   1  9000000 <int [5,443]>
#'   2  9500000 <int [59,183]>
#'   3 10000000 <int [45,040]>
#'   4 10500000 <int [96,504]>
#'
#' @return A list containing fragment lengths
#' @export
preprocess_frag_bed <- function(frag, bin_size) {
  chrom_list <- unique(seqnames(frag))
  result <- chrom_list %>%
    map(function(chrom) {
      frag <- frag[seqnames(frag) == chrom]
      fraglen_over_window(frag, bin_size = bin_size)
    })
  names(result) <- chrom_list
  result
}


#' @param block_size The genomic range to be processed will be divided into
#'   blocks, or segments, and assigned to workers. Each segment represents the
#'   smallest unit of the jobs. `block_size` controls the size of the segment.
#'   When there are more fragments, generally you need smaller `block_size` so
#'   that each job will not consume to much memory.
#' @export
contact_matrix <-
  function(fraglen_list,
           frag = NULL,
           bin_size = 500e3L,
           block_size = 10e6L,
           n_workers = 1L,
           subsample = 10e3L,
           min_sample_size = 100L,
           bootstrap = 1L,
           seed = NULL) {
    if (!is.null(frag)) {
      stop("Using frag as input is deprecated. Use fraglen_list instead.")
    }

    chrom <- names(fraglen_list)
    stopifnot(length(chrom) >= 1)

    result <- chrom %>% map(function(chrom) {
      fraglen <- fraglen_list[[chrom]]
      plan <- cofrag_plan(fraglen, bin_size = bin_size)
      calc_contact_matrix(
        fraglen,
        plan,
        n_workers = n_workers,
        block_size = block_size,
        stat_func = stat_func_ks(
          subsample = subsample,
          min_sample_size = min_sample_size,
          bootstrap = bootstrap
        ),
        bin_size = bin_size,
        seed = seed
      ) %>%
        build_gr(chrom = chrom, bin_size = bin_size)
    }) %>%
      do.call(c, args = .)
  }
