#!/usr/bin/env Rscript

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# # example
# script_args <- list(
#   input = here("sandbox/Pilot2-9.hg19.frag.bed.gz"),
#   output_dir = here("sandbox/"),
#   sample_id = "T1",
#   metrics = "ks",
#   genome = "hs37-1kg",
#   res = 500e3L,
#   ncores = 6L,
#   bootstrap = 3L,
#   subsample = 10e3L,
#   seed = 1228L,
#   chroms = c("20", "21", "22"),
#   exclude_chroms = c("18", "20"),
#   min_mapq = 30L,
#   min_fraglen = 50L,
#   max_fraglen = 500L,
#   intersect_region = NULL, # here("sandbox/cofrag/cfEW1.cna.neutral.bed"),
#   exclude_region = "encode.blacklist.hs37-1kg"
# )

if (interactive()) {
  library(tidyverse)
  library(magrittr)

  if (is.null(get0("script_args")))
    stop_quietly()
} else {
  # Run in CLI script mode
  parser <- optparse::OptionParser(
    option_list = list(
      optparse::make_option(c("-i", "--input")),
      optparse::make_option(c("-o", "--output-dir")),
      optparse::make_option(c("-s", "--sample-id")),
      optparse::make_option(c("--res"), type = "integer", default = 500e3L),
      optparse::make_option(c("--chroms"), default = NULL,
                            help = "Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used"),
      optparse::make_option(c("--method"), default = "NONE"),
      # optparse::make_option(
      #   c("--standard-hic"),
      #   type = "character",
      #   default = NULL,
      #   help = "A standard Hi-C map (.hic file) used for calculating HiCRep correlation scores"
      # ),
      # optparse::make_option(
      #   c("--standard-compartment"),
      #   default = NULL,
      #   help = "A standard compartment track (BED format) for calculating compartment correlation scores"
      # ),
      optparse::make_option(c("--juicer"), default = NULL,
                            help = "Path to the .jar file Juicer tools. If not provided, will download from Internet"),
      optparse::make_option(c("--java"), default = "java",
                            help = "Path to java")
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

  if (!is.null(script_args$method))
    script_args$method %<>% str_split(pattern = ":") %>% .[[1]]
  if (!is.null(script_args$chroms))
    script_args$chroms %<>% str_split(pattern = ":") %>% .[[1]]
  # if (!is.null(script_args$standard_hic))
  #   script_args$standard_hic %<>% str_split(pattern = ":") %>% .[[1]]
  # if (!is.null(script_args$standard_compartment))
  #   script_args$standard_compartment %<>% str_split(pattern = ":") %>% .[[1]]
}

# Build comment lines
comments <- c(
  paste0("cofragr version: ", as.character(packageVersion("cofragr"))),
  paste0("bedtorch version: ", as.character(packageVersion("bedtorch"))),
  # All items in script_args
  names(script_args) %>% purrr::map_chr(function(name) {
    v <- script_args[[name]]
    v_str <- if (!is.null(v))
      v %>% purrr::map_chr(as.character) %>% paste(collapse = ":")
    else
      ""
    paste0(name, "=", v_str)
  })
)

library(here)

logging::loginfo(str_interp("Argument summary:"))
comments %>% purrr::walk(function(v) logging::loginfo(v))






