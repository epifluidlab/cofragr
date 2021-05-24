#!/usr/bin/env Rscript

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


# # example
# script_args <- list(
#   input = here("sandbox/Pilot2-9.cofrag_cm.bed.gz"),
#   output_dir = here("sandbox/"),
#   sample_id = "T1",
#   res = 500e3L,
#   chroms = c("20", "21", "22"),
#   smooth = 3L,
#   method = "NONE",
#   juicer = NULL,
#   java = "java"
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
      # optparse::make_option(c("--method"), default = ""),
      # optparse::make_option(c("--smooth"), type = "integer", default = 3L, help = "Apply moving average on the compartment track. Sometimes this can remove some quirks and make the data more similar to Hi-C compartment tracks"),
      # optparse::make_option(
      #   c("--standard-hic"),
      #   type = "character",
      #   default = NULL,
      #   help = "A standard Hi-C map (.hic file) used for calculating HiCRep correlation scores"
      # ),
      optparse::make_option(
        c("--standard-compartment"),
        default = NULL,
        help = "A standard compartment track (BED format) for calculating compartment correlation scores"
      ),
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
  paste0("hictools version: ", as.character(packageVersion("hictools"))),
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

all_chroms <- system(str_interp("tabix -l ${script_args$input}"), intern = TRUE)
logging::loginfo(
  str_interp(
    "Found ${length(all_chroms)} chromosomes: ${paste(all_chroms, collapse = \",\")}"
  )
)

# Filter chroms
if (!is.null(script_args$chroms))
  all_chroms <- intersect(all_chroms, script_args$chroms)
if (!is.null(script_args$exclude_chroms))
  all_chroms <- setdiff(all_chroms, script_args$exclude_chroms)
logging::loginfo(
  str_interp(
    "Process ${length(all_chroms)} chromosomes: ${paste(all_chroms, collapse = \",\")}"
  )
)

comps <- expand_grid(
  chrom = all_chroms,
  smoothing = 1:3,
  method = c("juicer", "hictools")
) %>%
  pmap(function(chrom, smoothing, method) {
    hic_data <-
      hictools::load_hic_genbed(
        script_args$input,
        chrom = chrom,
        resol = script_args$res,
        type = "oe",
        norm = "NONE"
      )

    if (!is.null(script_args$standard_compartment))
      standard_comp <-
        bedtorch::read_bed(script_args$standard_compartment,
                           range = chrom,
                           use_gr = FALSE)
    else
      standard_comp <- "gene_density.hg19"

    if (method == "juicer") {
      comps <-
        hictools::compartment_juicer(hic_data, standard = standard_comp, smoothing = smoothing)
    } else {
      comps <-
        hictools::compartment_ht(hic_data,
                                 standard = standard_comp,
                                 type = "oe",
                                 smooth = smoothing)
    }

    comps <-
      cbind(comps,
            data.table::data.table(method = method, smooth = smoothing))
    comps[, score := ifelse(is.infinite(score) |
                              is.nan(score), NA, score)]
  }) %>%
  data.table::rbindlist(fill = TRUE)


system(str_interp("mkdir -p ${script_args$output_dir}"))
cofrag_comp_file <- str_interp("${script_args$output_dir}/${script_args$sample_id}.compartment.bed.gz")
logging::loginfo(str_interp("Write compartment scores: ${cofrag_comp_file}"))

bedtorch::write_bed(comps, file_path = cofrag_comp_file, comments = comments, tabix_index = FALSE)

logging::loginfo(str_interp("Analysis completed, with results located at ${script_args$output_dir}"))




