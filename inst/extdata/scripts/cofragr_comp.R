#!/usr/bin/env Rscript

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


# # example
# script_args <- list(
#   input = "~/Downloads/EE87929.cofrag_cm.bed.gz",
#   output_dir = here::here("sandbox/"),
#   sample_id = "EE87929",
#   res = 500e3L,
#   chroms = c("20", "21", "22"),
#   # method = c("juicer", "lieberman", "obs_exp"),
#   method = c("lieberman", "obs_exp"),
#   smooth = c(1L, 3L),
#   standard_compartment = "wbc",
#   genome = "hs37-1kg",
#   juicer = NULL,
#   java = "java"
# )


# script_args validity
validate_args <- function(args) {
  with(args, {
    # assert_that(
    #   is_null(standard_compartment) ||
    #     (
    #       is_scalar_character(standard_compartment) &&
    #         standard_compartment %in% c("wbc", "gene_density")
    #     )
    # )
    assert_that(is_scalar_integer(res) && res > 0)
    assert_that(is_scalar_character(genome) && genome %in% c("GRCh38", "hs37-1kg"))
    assert_that(length(method) >= 1 &&
                  all(method %in% c("juicer", "lieberman", "obs_exp")))
    # assert_that(length(smooth) >= 1 && is_integer(smooth) && all(smooth > 0))
  }) %>% invisible()
} 


if (interactive()) {
  library(tidyverse)
  library(magrittr)
  library(here)
  library(rlang)
  library(assertthat)

  if (is.null(get0("script_args")))
    stop_quietly()
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
      optparse::make_option(c("--res"), type = "integer", default = 500e3L),
      optparse::make_option(c("--chroms"), default = NULL,
                            help = "Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used"),
      optparse::make_option(c("--method"), default = "juicer:lieberman:obs_exp"),
      optparse::make_option(c("--smooth"), default = "1:3", help = "Apply moving average on the compartment track. Sometimes this can remove some quirks and make the data more similar to Hi-C compartment tracks"),
      optparse::make_option(
        c("--standard-compartment"),
        default = NULL,
        help = "A standard compartment track (BED format) for calculating compartment correlation scores. Can be either `wbc` or `gene_density`"
      ),
      optparse::make_option(
        c("-g", "--genome"),
        default = NULL,
        help = "Reference genome of the dataset. The default is hs37-1kg.",
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

  script_args$method %<>% str_split(pattern = ":") %>% .[[1]]
  script_args$smooth %<>% str_split(pattern = ":") %>% .[[1]] %>% map_int(as.integer)
  if (!is_null(script_args$chroms))
    script_args$chroms %<>% str_split(pattern = ":") %>% .[[1]]
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

# Standard compartment
if (is_null(script_args$standard_compartment)) {
  standard_comp <- NULL
} else{
  standard_comp <- local({
    res_kbp <- as.integer(script_args$res / 1e3L)
    comp_name <- script_args$standard_compartment
    if (comp_name == "gene_density")
      comp_name <-
      str_interp("gencode.v30.b37.gene_density.${res_kbp}kbp")
    else if (comp_name == "wbc")
      comp_name <-
      str_interp("wbc.rep1.compartment.hs37-1kg.NONE.${res_kbp}kbp")
    
    logging::loginfo(str_interp("Loading standard compartment: ${comp_name}"))
    env <- rlang::env()
    data(list = comp_name, envir = env, package = "hictools")
    env[[comp_name]]
  })
}

cofrag_comp_grid <- expand_grid(
  chrom = all_chroms,
  method = script_args$method,
  smooth = script_args$smooth,
  # metric = c("score", "hellinger")
  metric = "score"
) 

cofrag_comp_results <- cofrag_comp_grid %>%
  pmap(function(chrom, method, smooth, metric) {
    logging::loginfo(str_interp("Processing: chrom:${chrom} ${method} smooth:${smooth} ${metric}"))
    hic_data <-
      hictools::load_hic_genbed(
        script_args$input,
        chrom,
        resol = script_args$res,
        type = "observed",
        norm = "NONE",
        genome = script_args$genome,
        score_col = metric,
        scale_score = TRUE
      )
    
    comps <-
      hictools::get_compartment(hic_data,
                                method = method,
                                standard = standard_comp,
                                smooth = smooth,
                                genome = script_args$genome)
    
    comps$method <- method
    comps$smooth <- smooth
    comps$metric <- metric
    
    comps$score <- local({
      score <- comps$score
      ifelse(is.infinite(score) | is.na(score), NA, score)
    })
    
    seqinfo <- bedtorch::get_seqinfo(script_args$genome)
    suppressWarnings({
      GenomeInfoDb::seqlevels(comps) <- GenomeInfoDb::seqlevels(seqinfo)
      GenomeInfoDb::seqinfo(comps) <- bedtorch::get_seqinfo(script_args$genome)
      GenomicRanges::trim(comps)
    })
    
    if (!is.null(standard_comp))
      correlation <- hictools::comp_correlation(comps, standard_comp)
    else
      correlation <- NA

    list(comps = comps, correlation = correlation)
  })

comps <- cofrag_comp_results %>% map(~ .$comps) %>% do.call(what = c, args = .)

# cofrag_correlation <-
#   map2_dfr(seq.int(nrow(cofrag_comp_grid)), 
#            cofrag_comp_results, 
#            function(grid_idx, result) {
#     c(cofrag_comp_grid[grid_idx,], result$correlation)
#   })

system(str_interp("mkdir -p ${script_args$output_dir}"))
cofrag_comp_file <- str_interp("${script_args$output_dir}/${script_args$sample_id}.compartment.bed.gz")
cofrag_correlation_file <- str_interp("${script_args$output_dir}/${script_args$sample_id}.compartment_correlation.tsv")
logging::loginfo(str_interp("Write compartment scores: ${cofrag_comp_file}"))
logging::loginfo(str_interp("Write compartment correlation: ${cofrag_correlation_file}"))

bedtorch::write_bed(comps, file_path = cofrag_comp_file, comments = comments, tabix_index = FALSE)
# write_tsv(cofrag_correlation, file = cofrag_correlation_file)

logging::loginfo(str_interp("Analysis completed, with results located at ${script_args$output_dir}"))




