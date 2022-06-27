#' Read a fragment data file
#'
#' A fragment data file is essentially just a BED file. The reason we need a
#' special function here to load a fragment data file is that, by convention,
#' the data frame should contain a column representing MAPQ scores. However, the
#' data file may not have a column header. Therefore, we may need to guess which
#' column represents MAPQ scores.
#' @param file_path a character vector representing paths of fragment files.
#'   This argument accepts multiple input files which is convenient for sample
#'   pooling.
#' @param range read the entire file if `range` is `NULL`.
#' @export
read_fragments <- function(file_path, range = NULL, genome = NULL) {
  frag <- lapply(file_path, function(path) {
    bedtorch::read_bed(file_path = path,
                       range = range,
                       genome = genome)
  }) %>%
    do.call(what = c, args = .)
  GenomicRanges::strand(frag) <- "*"
  frag_metadata <- mcols(frag)
  
  if ("mapq" %in% colnames(frag_metadata)) {
    if (!is.integer(frag_metadata$mapq)) {
      logging::logerror("MAPQ should be integers")
      print(frag)
      stop()
    }
  } else if (ncol(frag_metadata) >= 2 && is.integer(frag_metadata[[2]])) {
    frag_metadata$mapq <- frag_metadata[[2]]
  } else {
    logging::logwarn("Cannot find MAPQ in the data frame")
  }
  
  # Only need mapq. Drop other columns
  mcols(frag) <-
    frag_metadata[intersect(c("mapq", "cigar1", "cigar2"), colnames(frag_metadata))]
  
  if ("cigar1" %in% colnames(frag_metadata)) {
    cigar_select <-
      grepl(pattern = "^[0-9]+M", frag$cigar1) &
      grepl(pattern = "[0-9]+M$", frag$cigar2)
    frag <- frag[cigar_select]
  }
  
  # Each fragment is represented using the mid-point
  frag$len <- GenomicRanges::width(frag)
  mid_point <- as.integer(round((GenomicRanges::start(frag) + GenomicRanges::end(frag)) / 2))
  GenomicRanges::start(frag) <- mid_point
  GenomicRanges::end(frag) <- mid_point
  frag <- sort(frag)
  
  print(frag)
  
  return(frag)
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
  stopifnot(bin_size %in% c(50e3L, 100e3L, 250e3L, 500e3L, 1e6L, 2.5e6L, 5e6L))

  # Only allow single chromosome fragment data
  chrom <- unique(seqnames(frag))
  stopifnot(length(chrom) == 1)

  # # Find mid-points
  # midpoints <-
  #   as.integer(round(start(frag) + (width(frag) - 1) / 2))
  # frag$len <- width(frag)
  # start(frag) <- midpoints
  # width(frag) <- 1

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
    if (isTRUE(min(n_frag1, n_frag2) < min_sample_size))
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
    # p_value <- mean(v)
    # p_value_sd <- sd(v)
    
    
    # Calculate hellinger distance
    if (identical(len1, len2)) {
      score_hellinger <- 1
    } else {
      bkde1 <- pmax(philentropy::binned.kernel.est(data = len1, range.data = c(1, 1000), gridsize = 1000)$y, 0)
      bkde2 <- pmax(philentropy::binned.kernel.est(data = len2, range.data = c(1, 1000), gridsize = 1000)$y, 0)
      hellinger_dist <- philentropy::distance(
        x = rbind(bkde1, bkde2),
        method = "hellinger",
        est.prob = "empirical",
        mute.message = TRUE
      ) / 2
      score_hellinger <- 1 - hellinger_dist
    }

    # score, n_frag1 and n_frag2 doesn't change across different bootstrap iterations
    # To save space, we only record these values at bootstrap #1. For other iterations,
    # we write NA instead.
    # For p-values, however, we record all of them faithfully since they vary across iterations.
    list(
      score = c(score, rep(NA, bootstrap - 1)),
      n_frag1 = c(n_frag1, rep(NA, bootstrap - 1)),
      n_frag2 = c(n_frag2, rep(NA, bootstrap - 1)),
      bootstrap = seq.int(bootstrap),
      p_value = v,
      hellinger = c(score_hellinger, rep(NA, bootstrap - 1))
    )
  }
}


calc_contact_matrix <- function(fraglen,
                                cofrag_plan,
                                stat_func,
                                bin_size = 500e3L,
                                block_size = 10e6L,
                                n_workers = 1L,
                                parallel = FALSE,
                                seed = NULL) {
  # If not in parallel mode, need just one core
  if (!parallel)
    n_workers <- 1

  if (parallel && n_workers > 1) {
    # cl <- makeForkCluster(n_workers, outfile = "")
    cl <- makeCluster(n_workers, outfile = "")
    registerDoParallel(cl)
    logging::loginfo(str_interp("Successfully registered a cluster with ${n_workers} workers"))
    on.exit(stopCluster(cl), add = TRUE)
    
    `%loop%` <- doRNG::`%dorng%`
  } else
    `%loop%` <- foreach::`%do%`

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
      fraglen_map <- list()
      unique(c(id_list_1, id_list_2)) %>%
        walk(function(id) {
          fraglen_map[[as.character(id)]] <<- fraglen$len[[id]]
        })

      list(
        block_id = idx,
        id_list_1 = id_list_1,
        id_list_2 = id_list_2,
        fraglen_map = fraglen_map
      )
    })

  gc()
  logging::loginfo(str_interp("Completed building plans. Current memory usage: ${as.numeric(lobstr::mem_used()) / 1e6}"))

  if (!is.null(seed))
    set.seed(seed)
  
  foreach(
    data = foreach_layout,
    .export = "n_blocks",
    .combine = "rbind",
    .multicombine = TRUE
  ) %loop% {
    block_id <- data$block_id

    id_list_1 <- data$id_list_1
    id_list_2 <- data$id_list_2
    fraglen_map <- data$fraglen_map

    worker_name <-
      paste(Sys.info()[['nodename']], Sys.getpid(), sep = '-')
    mem_mb <- as.numeric(lobstr::mem_used()) / 1e6
    logging::loginfo(stringr::str_interp("${worker_name}: processing ${block_id}/${n_blocks}, memory used: ${mem_mb} MB"))


    stopifnot(length(id_list_1) == length(id_list_2))

    purrr::map_dfr(1:length(id_list_1), function(idx) {
      id.1 <- id_list_1[idx]
      id.2 <- id_list_2[idx]

      len1 <- fraglen_map[[as.character(id.1)]]
      len2 <- fraglen_map[[as.character(id.2)]]

      result <- stat_func(len1, len2)
      if (!is.null(result)) {
        result$id.1 <- id.1
        result$id.2 <- id.2
      }
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
           parallel = FALSE,
           subsample = 10e3L,
           min_sample_size = 100L,
           bootstrap = 1L,
           seed = NULL) {
    assert_that(is_list(fraglen_list))
    assert_that(is_null(frag), msg = "Using frag as input is deprecated. Use fraglen_list instead.")
    bin_size <- as.integer(bin_size)
    assert_that(is_scalar_integer(bin_size) && bin_size %in% c(50e3L, 100e3L, 250e3L, 500e3L, 1e6L, 2.5e6L, 5e6L), msg = "Invalid bin_size")
    block_size <- as.integer(block_size)
    assert_that(is_scalar_integer(block_size) && block_size > 0 && block_size %% bin_size == 0)
    assert_that(is_scalar_integer(n_workers) && n_workers >= 1)
    assert_that(is_scalar_integer(subsample) && subsample >= 100)
    assert_that(is_null(min_sample_size) || (is_scalar_integer(min_sample_size) && min_sample_size > 0))
    assert_that(is_scalar_integer(bootstrap) && bootstrap >= 1)
    assert_that(is_null(seed) || is_scalar_integer(seed))
    
    if (!is.null(frag)) {
      stop("Using frag as input is deprecated. Use fraglen_list instead.")
    }

    chrom <- names(fraglen_list)
    stopifnot(length(chrom) >= 1)

    chrom %>% map(function(chrom) {
      fraglen <- fraglen_list[[chrom]]
      plan <- cofrag_plan(fraglen, bin_size = bin_size)
      cm <- calc_contact_matrix(
        fraglen,
        plan,
        n_workers = n_workers,
        parallel = parallel,
        block_size = block_size,
        stat_func = stat_func_ks(
          subsample = subsample,
          min_sample_size = min_sample_size,
          bootstrap = bootstrap
        ),
        bin_size = bin_size,
        seed = seed
      )
      cm %>% rename(
        pos1 = start1,
        pos2 = start2
      ) %>% mutate(
        chrom1 = chrom,
        chrom2 = chrom
      ) %>%
        select(chrom1, pos1, chrom2, pos2, everything()) %>%
        hictools::ht_table(resol = bin_size, type = "cofrag", norm = "NONE")
    }) %>%
      hictools::concat_hic()
  }
