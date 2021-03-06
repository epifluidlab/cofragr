#!/usr/bin/env Rscript

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# script_args validity
validate_args <- function(args) {
  with(args, {
    assert_that(
      is_scalar_character(metrics) &&
        metrics %in% c("ks", "hellinger"),
      msg = paste0("Invalid metrics: ", metrics)
    )
    assert_that(
      is_scalar_character(genome) &&
        genome %in% c("hs37-1kg", "GRCh37", "GRCh38"),
      msg = str_interp("Unsupported genome: ${genome}")
    )
    assert_that(is_scalar_integer(res) && res > 0)
    assert_that(is_scalar_integer(block_size) &&
      block_size > 0 &&
      block_size %% res == 0,
    msg = "block_size must be a positive integer which is the multiple of res"
    )
    assert_that(is_scalar_integer(ncores) && ncores >= 1)
    assert_that(is_scalar_integer(bootstrap) && bootstrap >= 1)
    assert_that(is_scalar_integer(subsample) &&
      subsample >= 100, msg = "subsample must be a positive integer no less than 100")
    assert_that(is_scalar_integer(min_mapq) && min_mapq >= 0)
    assert_that(is_scalar_integer(min_fraglen) && min_fraglen >= 1)
    assert_that(is_scalar_integer(max_fraglen) && max_fraglen >= 1)
    assert_that(is_null(seed) || is_scalar_integer(seed))
    assert_that(!isTRUE(max_fraglen <= min_fraglen), msg = "max_fraglen must be larger than min_fraglen")
  }) %>% invisible()
}


# example
script_args <- list(
  input = "~/dev/gateway_project/sandbox/frag/Pilot2_9.hg19.frag.bed.gz",
  output_dir = here::here("sandbox/"),
  sample_id = "Pilot2_9.hg19",
  metrics = "ks",
  genome = "hs37-1kg",
  res = 500e3L,
  block_size = 10e6L,
  ncores = 2L,
  bootstrap = 5L,
  subsample = 10e3L,
  seed = 1228L,
  chroms = c("21"),
  exclude_chroms = NULL,
  min_mapq = 30L,
  min_fraglen = 100L,
  max_fraglen = 350L,
  intersect_region = NULL, # here("sandbox/cofrag/cfEW1.cna.neutral.bed"),
  exclude_region = "encode.blacklist",
  disable_parallel = FALSE
)


if (interactive()) {
  library(tidyverse)
  library(magrittr)
  library(here)
  library(rlang)
  library(assertthat)

  if (is.null(get0("script_args"))) {
    stop_quietly()
  }

  # Build a string to record all parameters
  param_str <-
    paste(names(script_args),
      script_args,
      sep = "=",
      collapse = "; "
    )
} else {
  library(here)
  library(rlang)
  library(assertthat)

  # Run in CLI script mode
  parser <- optparse::OptionParser(
    option_list = list(
      optparse::make_option(c("-i", "--input")),
      optparse::make_option(c("-o", "--output-dir")),
      optparse::make_option(c("-s", "--sample-id")),
      optparse::make_option(c("-m", "--metrics"),
        default = "ks",
        help = "Metrics for calculating contact scores. Must be ks or hellinger [ks]"
      ),
      optparse::make_option(c("-g", "--genome"),
        default = "GRCh38",
        help = "Reference genome for the input. Must be hs37-1kg, GRCh37, or GRCh38 [GRCh38]"
      ),
      optparse::make_option(c("--res"),
        type = "integer", default = 500e3L,
        help = "Resolution, aka bin size [500000]"
      ),
      optparse::make_option(c("--standard-compartment"),
        type = "character", default = NULL,
        help = "Path to a BED file for standard compartment scores [NULL]"
      ),
      optparse::make_option(c("--block-size"),
        type = "integer", default = 10e6L,
        help = "Block size. Must be multiple of --res [10000000]"
      ),
      optparse::make_option(c("-n", "--ncores"),
        type = "integer", default = 1L,
        help = "Number of threads [1]"
      ),
      optparse::make_option(c("--bootstrap"),
        type = "integer", default = 1L,
        help = "Number of bootstrap iterations [1]"
      ),
      optparse::make_option(c("--subsample"),
        type = "integer", default = NULL,
        help = "Subsample size. NULL means do not perform subsampling [NULL]"
      ),
      optparse::make_option(c("--seed"),
        type = "integer", default = NULL,
        help = "Seed for random number generator [NULL]"
      ),
      optparse::make_option(c("--chroms"),
        default = NULL,
        help = "Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used [NULL]"
      ),
      optparse::make_option(c("--exclude-chroms"),
        default = NULL,
        help = "Exclude chromosomes from the analysis. Separated by colons, such as 12:16:X [NULL]"
      ),
      optparse::make_option(
        c("--min-mapq"),
        type = "integer",
        default = 1L,
        help = "Minimal MAPQ for fragments to take part in the analysis [1]"
      ),
      optparse::make_option(
        c("--min-fraglen"),
        type = "integer",
        default = 100L,
        help = "Minimal length for fragments to take part in the analysis [100]"
      ),
      optparse::make_option(
        c("--max-fraglen"),
        type = "integer",
        default = 350L,
        help = "Maximal length for fragments to take part in the analysis [300]"
      ),
      optparse::make_option(
        c("--intersect-region"),
        type = "character",
        default = NULL,
        help = "BED files defining regions to be intersected with the fragment data file, separated by colon [NULL]"
      ),
      optparse::make_option(
        c("--exclude-region"),
        type = "character",
        default = NULL,
        # default = "encode.blacklist",
        help = "BED files defining regions to be excluded from the analysis, separated by colon [NULL]"
      ),
      optparse::make_option(c("--disable-parallel"), action = "store_true", default = FALSE, help = "Run in non-parallel mode [FALSE]")
    )
  )
  script_args <-
    optparse::parse_args(
      parser,
      args = commandArgs(trailingOnly = TRUE),
      convert_hyphens_to_underscores = TRUE
    )

  library(tidyverse)
  library(magrittr)

  if (is.null(script_args$seed)) {
    script_args$seed <- as.integer(round(runif(1) * 1e6L))
    logging::loginfo(str_interp("Randomly pick a seed: ${script_args$seed}"))
  }
  if (!is.null(script_args$chroms)) {
    script_args$chroms %<>% str_split(pattern = ":") %>% .[[1]]
  }
  if (!is.null(script_args$exclude_chroms)) {
    script_args$exclude_chroms %<>% str_split(pattern = ":") %>% .[[1]]
  }
  if (!is.null(script_args$intersect_region)) {
    script_args$intersect_region %<>% str_split(pattern = ":") %>% .[[1]]
  }
  if (!is.null(script_args$exclude_region)) {
    script_args$exclude_region %<>% str_split(pattern = ":") %>% .[[1]]
  }

  # Build a string to record all parameters
  param_str <- paste0("Parameters: ", paste(commandArgs(trailingOnly = TRUE), collapse = " "))
}

validate_args(script_args)

# Build comment lines
comments <- c(
  paste0("cofragr_version=", as.character(packageVersion("cofragr"))),
  paste0("hictools_version=", as.character(packageVersion("hictools"))),
  paste0("bedtorch_version=", as.character(packageVersion("bedtorch"))),
  # All items in script_args
  names(script_args) %>% purrr::map_chr(function(name) {
    v <- script_args[[name]]
    v_str <- if (!is.null(v)) {
      v %>%
        purrr::map_chr(as.character) %>%
        paste(collapse = ":")
    } else {
      ""
    }
    paste0(name, "=", v_str)
  })
)

logging::loginfo(str_interp("Argument summary:"))
comments %>% purrr::walk(function(v) logging::loginfo(v))

all_chroms <- system(str_interp("tabix -l ${script_args$input}"), intern = TRUE)
logging::loginfo(
  str_interp(
    "Found ${length(all_chroms)} chromosomes: ${paste(all_chroms, collapse = \",\")}"
  )
)
# Filter chroms
if (!is_null(script_args$chroms)) {
  all_chroms <- intersect(all_chroms, script_args$chroms)
}
if (!is_null(script_args$exclude_chroms)) {
  all_chroms <- setdiff(all_chroms, script_args$exclude_chroms)
}
logging::loginfo(
  str_interp(
    "Process ${length(all_chroms)} chromosomes: ${paste(all_chroms, collapse = \",\")}"
  )
)


# Users may provide multiple files for excluded regions, and for each file the regions may be overlapping.
# Merge all excluded regions together.
merge_regions <- function(regions, genome = NULL) {
  if (is_null(regions)) {
    return(NULL)
  }

  regions %>%
    map(function(r) {
      # Use regions shipped with the package
      if (r == "encode.blacklist") {
        local({
          env <- rlang::env()
          dataset_name <- paste0("wgEncodeDacMapabilityConsensusExcludable.", genome)
          data(list = dataset_name, package = "cofragr", envir = env)
          region <- env[[dataset_name]]
          assert_that(!is_null(region))
          region
        })
      } else {
        bedtorch::read_bed(r, genome = genome)
      }
    }) %>%
    do.call(what = c, args = .) %>%
    bedtorch::merge_bed()
}


filter_mapq <- function(frag, min_mapq) {
  if (length(frag) == 0) {
    logging::loginfo("filter_mapq: empty fragments, skipping")
    return(frag)
  }

  if (!is.null(min_mapq)) {
    frag <- frag[frag$mapq >= min_mapq]
    # logging::loginfo(str_interp("# of fragments after MAPQ filter: ${length(frag)}"))
  }
  frag
}

filter_fraglen <- function(frag, min_fraglen = NULL, max_fraglen = NULL) {
  if (length(frag) == 0) {
    logging::loginfo("filter_length: empty fragments, skipping")
    return(frag)
  }

  if (is.null(min_fraglen) && is.null(max_fraglen)) {
    return(frag)
  } else if (is.null(min_fraglen)) {
    frag <- frag[frag$len <= max_fraglen]
  } else if (is.null(max_fraglen)) {
    frag <- frag[frag$len >= min_fraglen]
  } else {
    frag <- frag[between(frag$len, min_fraglen, max_fraglen)]
  }

  frag
}

cofrag_cm <- all_chroms %>%
  map(function(chrom) {
    logging::loginfo(str_interp("Loading fragments for chromosome ${chrom}"))

    frag <- cofragr::read_fragments(file_path = script_args$input, range = chrom, genome = script_args$genome)

    # if ("cigar1" %in% frags_cols && "cigar2" %in% frags_cols) {
    #   # Exclude soft clipping
    #   frags %<>%
    #     filter(!(str_detect(cigar1, pattern = "^[0-9]+S") | str_detect(cigar2, pattern = "[0-9]+S$")))
    # }

    logging::loginfo(str_interp("# of fragments loaded for chr${chrom}: ${length(frag)}"))

    if (!is_null(script_args$intersect_region)) {
      intersect_region <- merge_regions(script_args$intersect_region, genome = script_args$genome)
      frag %<>% bedtorch::intersect_bed(intersect_region, mode = "unique")
      logging::loginfo(str_interp("# of fragments after intersecting regions: ${length(frag)}"))
    }

    if (!is_null(script_args$exclude_region)) {
      exclude_region <- merge_regions(script_args$exclude_region, genome = script_args$genome)
      frag %<>% bedtorch::exclude_bed(exclude_region)
      logging::loginfo(str_interp("# of fragments after excluding regions: ${length(frag)}"))
    }

    if (!is_null(script_args$min_mapq)) {
      frag %<>% filter_mapq(min_mapq = script_args$min_mapq)
      logging::loginfo(str_interp("# of fragments after filtering by MAPQ: ${length(frag)}"))
    }

    if (!is_null(script_args$min_fraglen) || !is_null(script_args$max_fraglen)) {
      frag %<>% filter_fraglen(min_fraglen = script_args$min_fraglen, max_fraglen = script_args$max_fraglen)
      logging::loginfo(str_interp("# of fragments after filtering by fragment length: ${length(frag)}"))
    }

    # Minimal # of fragments: 10,000
    if (length(frag) < 10e3L) {
      logging::logwarn(str_interp("Not enough fragments to call the matrix: ${length(frag)}"))
      return(NULL)
    }

    fraglen_list <- cofragr::preprocess_frag_bed(frag, bin_size = script_args$res)
    # Explicitly remove frag to save memory
    rm(frag)
    gc()

    logging::loginfo(str_interp("Completed preprocessing. Current memory usage: ${as.numeric(lobstr::mem_used()) / 1e6}"))

    cofragr::contact_matrix(
      fraglen_list,
      bin_size = script_args$res,
      block_size = script_args$block_size,
      n_workers = script_args$ncores,
      parallel = !script_args$disable_parallel,
      subsample = script_args$subsample,
      min_sample_size = 100L,
      bootstrap = script_args$bootstrap,
      seed = script_args$seed
    )
  }) %>%
  hictools::concat_hic()


output_dir <- script_args$output_dir
system(str_interp("mkdir -p ${output_dir}"))
cofrag_cm_file <- str_interp("${output_dir}/${script_args$sample_id}.cofrag_cm.bed.gz")
logging::loginfo(str_interp("Write cofragmentation matrix: ${cofrag_cm_file}"))

hictools::write_hic_bedtorch(cofrag_cm, file_path = cofrag_cm_file, comments = comments)
# cofragr::write_contact_matrix(cofrag_cm, file_path = cofrag_cm_file, comments = comments)

logging::loginfo(str_interp("Analysis completed, with results located at ${output_dir}"))