# MIT License
#
# Copyright (c) 2021 Haizi Zheng
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Author: Haizi Zheng
# Copyright: Copyright 2021, Haizi Zheng
# Email: haizi.zh@gmail.com
# License: MIT
#
# Perform compartment analysis

library(here)
library(rlang)
library(tidyverse)
library(magrittr)

script_version <- "0.1.2"

invisible(
  assertthat::assert_that(
    packageVersion("bedtorch") >= "0.1.12.11",
    packageVersion("hictools") >= "0.1.16"
  )
)

if (interactive()) {
  input_file <- here("sandbox/EE86274.cofrag_cm.bed.gz")
  script_options <- list(
    output = here("sandbox/EE86274.compartment.fanc.bed.gz"),
    res = 500e3L,
    chroms = as.character(17:22),
    method = "juicer",
    smooth = c(1L, 3L),
    reference = here("cofrag/data/wbc_rep1/wbc_rep1.juicer_comps.KR.500kbp.bed"),
    genome = "GRCh37",
    num_cores = 6L
  )
} else {
  # Run in CLI script mode
  parser <- optparse::OptionParser(
    option_list = list(
      optparse::make_option(
        c("-o", "--output", type = "character", help = "Output file")
      ),
      optparse::make_option(
        c("--res"),
        type = "integer",
        default = 500e3L,
        help = "Resolution used in the compartment analysis [500000]"
      ),
      optparse::make_option(c("--chroms"),
        default = NULL,
        help = "Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used [NULL]"
      ),
      optparse::make_option(c("--method"), default = "juicer", help = "Methods for compartment analysis. Must be juicer or fanc [juicer]"),
      optparse::make_option(c("--smooth"), default = "1,3", help = "Apply moving average on the compartment track. Sometimes this can remove some quirks and make the data more similar to Hi-C compartment tracks [1,3]"),
      optparse::make_option(c("--reference"),
        default = NULL,
        help = "A reference tandard compartment track (BED format) for calculating compartment correlation scores. Can be `gc`, `gene.density` or a file path [NULL]"
      ),
      optparse::make_option(c("-g", "--genome"),
        default = NULL,
        help = "Reference genome of the dataset [NULL]",
      ),
      optparse::make_option(c("-n", "--num-cores"), default = 1L, help = "Number of cores for parallel computation [1]")
    )
  )

  script_args <-
    optparse::parse_args(
      parser,
      args = commandArgs(trailingOnly = TRUE),
      convert_hyphens_to_underscores = TRUE,
      positional_arguments = TRUE
    )

  input_file <- script_args$args
  invisible(assertthat::assert_that(rlang::is_scalar_character(input_file)))
  script_options <- script_args$options

  script_options$smooth %<>% str_split(pattern = ",") %>%
    .[[1]] %>%
    map_int(as.integer)
  if (!is_null(script_options$chroms)) {
    script_options$chroms %<>% str_split(pattern = ",") %>% .[[1]]
  }
  script_options$res <- as.integer(script_options$res)
  script_options$num_cores <- as.integer(script_options$num_cores)
}

# Build comment lines
comments <- c(
  paste0("script_version=", script_version),
  paste0("timestamp=", lubridate::now() %>% format("%Y-%m-%dT%H:%M:%S%z")),
  # All items in options
  setdiff(names(script_options), "help") %>%
    purrr::map_chr(function(name) {
      v <- script_options[[name]]
      v_str <- if (!is.null(v)) {
        v %>%
          purrr::map_chr(as.character) %>%
          paste(collapse = ",")
      } else {
        ""
      }
      paste0(name, "=", v_str)
    })
)

logging::loginfo(str_interp("Argument summary:"))
comments %>% purrr::walk(function(v) logging::loginfo(v))

all_chroms <- system(str_interp("tabix -l ${input_file}"), intern = TRUE)
logging::loginfo(
  str_interp(
    "Found ${length(all_chroms)} chromosomes: ${paste(all_chroms, collapse = \",\")}"
  )
)

# Filter chroms
if (!is.null(script_options$chroms)) {
  all_chroms <- intersect(all_chroms, script_options$chroms)
}
if (!is.null(script_options$exclude_chroms)) {
  all_chroms <- setdiff(all_chroms, script_options$exclude_chroms)
}
logging::loginfo(
  str_interp(
    "Process ${length(all_chroms)} chromosomes: ${paste(all_chroms, collapse = \",\")}"
  )
)

if (script_options$reference == "gc") {
  ref_comp <- local({
    env <- new.env()
    res <- paste0(script_options$res %/% 1000, "kbp")
    track_name <- stringr::str_interp("gc.${script_options$genome}.${res}")
    data(
      list = track_name,
      package = "hictools",
      envir = env
    )
    env[[track_name]]
  })
} else if (assertthat::is.readable(script_options$reference)) {
  ref_comp <-
    bedtorch::read_bed(
      file_path = script_options$reference,
      genome = script_options$genome
    )
} else if (!is_null(script_options$reference)) {
  ref_comp <- local({
    env <- new.env()
    data(
      list = script_options$reference,
      package = "hictools",
      envir = env
    )
    env[[track_name]]
  })
} else {
  ref_comp <- NULL
}


if (script_options$method == "juicer") {
  if (script_options$num_cores == 1) {
    bpparam <- BiocParallel::bpparam("SerialParam")
  } else {
    bpparam <- BiocParallel::MulticoreParam(workers = script_options$num_cores)
  }

  comps <-
    BiocParallel::bplapply(
      all_chroms,
      BPPARAM = bpparam,
      FUN = function(chrom) {
        hic_data <-
          hictools::load_hic(
            input_file,
            chrom = chrom,
            resol = script_options$res,
            genome = script_options$genome,
            type = "cofrag"
          )

        hictools::get_compartment(
          hic_matrix = hic_data,
          chrom = chrom,
          method = script_options$method,
          reference = ref_comp
        )
      }
    ) %>%
    do.call(what = c, args = .)

  comps$smooth <- 1L
  comps_smoothed <- setdiff(script_options$smooth, 1) %>%
    map(function(smooth) {
      unique(GenomicRanges::seqnames(comps)) %>%
        map(function(chrom) {
          comps <- comps[GenomicRanges::seqnames(comps) == chrom]
          comps$smooth <- smooth

          comps$score <- zoo::rollmean(
            comps$score,
            k = smooth,
            na.pad = TRUE,
            na.rm = TRUE,
            align = "center"
          )
          comps
        }) %>%
        do.call(what = c, args = .)
    }) %>%
    do.call(what = c, args = .)
  comps <- c(comps, comps_smoothed)

  comps$score <- local({
    score <- comps$score
    ifelse(is.infinite(score) | is.na(score), NA, score)
  })
  comps <- comps[!is.na(comps$score)]
} else {
  hic_data <-
    hictools::load_hic(
      input_file,
      chrom = all_chroms,
      resol = script_options$res,
      genome = script_options$genome
    )

  comps <- hictools::get_compartment(
    hic_matrix = hic_data,
    method = script_options$method,
    reference = ref_comp
  )

  comps$smooth <- 1L
  comps_smoothed <- setdiff(script_options$smooth, 1) %>%
    map(function(smooth) {
      unique(GenomicRanges::seqnames(comps)) %>%
        map(function(chrom) {
          comps <- comps[GenomicRanges::seqnames(comps) == chrom]
          comps$smooth <- smooth

          comps$PC1 <- zoo::rollmean(
            comps$PC1,
            k = smooth,
            na.pad = TRUE,
            na.rm = TRUE,
            align = "center"
          )

          comps$PC2 <- zoo::rollmean(
            comps$PC2,
            k = smooth,
            na.pad = TRUE,
            na.rm = TRUE,
            align = "center"
          )
          comps
        }) %>%
        do.call(what = c, args = .)
    }) %>%
    do.call(what = c, args = .)
  comps <- c(comps, comps_smoothed)

  comps$PC1 <- local({
    PC1 <- comps$PC1
    ifelse(is.infinite(PC1) | is.na(PC1), NA, PC1)
  })
  comps$PC2 <- local({
    PC2 <- comps$PC2
    ifelse(is.infinite(PC2) | is.na(PC2), NA, PC2)
  })
  comps <- comps[!is.na(comps$PC1) | !is.na(comps$PC2)]
}

comps$method <- script_options$method
comps <- bedtorch::as.bedtorch_table(comps)
data.table::setorder(comps, "chrom", "start")
# comps %<>% arrange(chrom, start)
bedtorch::write_bed(comps, file_path = script_options$output, comments = comments)