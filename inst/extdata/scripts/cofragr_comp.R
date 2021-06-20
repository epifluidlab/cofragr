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


# script_args validity
validate_args <- function(args) {
  with(args, {
    assert_that(
      is_null(standard_compartment) ||
        (
          is_scalar_character(standard_compartment) &&
            standard_compartment %in% c("wbc", "gene_density")
        )
    )
    assert_that(is_scalar_integer(res) && res > 0)
    assert_that(is_null(genome) ||
                  (is_scalar_character(genome) && genome == "hs37-1kg"))
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
      optparse::make_option(c("--method"), default = "juicer:lieberman"),
      # optparse::make_option(c("--smooth"), default = "1:3", help = "Apply moving average on the compartment track. Sometimes this can remove some quirks and make the data more similar to Hi-C compartment tracks"),
      # optparse::make_option(
      #   c("--standard-hic"),
      #   type = "character",
      #   default = NULL,
      #   help = "A standard Hi-C map (.hic file) used for calculating HiCRep correlation scores"
      # ),
      optparse::make_option(
        c("--standard-compartment"),
        default = NULL,
        help = "A standard compartment track (BED format) for calculating compartment correlation scores. Can be wbc or gene_density"
      ),
      optparse::make_option(
        c("-g", "--genome"),
        default = "hs37-1kg",
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
  # script_args$smooth %<>% str_split(pattern = ":") %>% .[[1]] %>% map_int(as.integer)
  if (!is_null(script_args$chroms))
    script_args$chroms %<>% str_split(pattern = ":") %>% .[[1]]
}

validate_args(script_args)

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

comps <- expand_grid(
  chrom = all_chroms,
  method = script_args$method
) %>%
  pmap(function(chrom, method) {
    hic_data <-
      hictools::load_hic_genbed(
        script_args$input,
        chrom = chrom,
        resol = script_args$res,
        type = "cofrag",
        norm = "NONE",
        genome = script_args$genome
      )
    comps <-
      hictools::get_compartment(hic_data,
                                method = method,
                                standard = standard_comp,
                                smoothing = 1,
                                genome = script_args$genome)
    comps$method <- method
    
    comps <- 1:3 %>%
      map(function(smoothing) {
        comps$smooth <- smoothing
        if (smoothing > 1)
          comps$score <- zoo::rollmean(comps$score, k = smoothing, na.pad = TRUE, na.rm = TRUE, align = "center")
        comps[!is.na(comps$score)]
      }) %>%
      do.call(c, args = .)
    
    # comps$score <- local({
    #   score <- comps$score
    #   ifelse(is.infinite(score) | is.na(score), NA, score)
    # })
    
    seqinfo <- bedtorch::get_seqinfo(script_args$genome)
    suppressWarnings({
      GenomeInfoDb::seqlevels(comps) <- GenomeInfoDb::seqlevels(seqinfo)
      GenomeInfoDb::seqinfo(comps) <- bedtorch::get_seqinfo(script_args$genome)
      GenomicRanges::trim(comps)
    })
  })

comps <- rlang::exec(c, !!!comps)

system(str_interp("mkdir -p ${script_args$output_dir}"))
cofrag_comp_file <- str_interp("${script_args$output_dir}/${script_args$sample_id}.compartment.bed.gz")
logging::loginfo(str_interp("Write compartment scores: ${cofrag_comp_file}"))

bedtorch::write_bed(comps, file_path = cofrag_comp_file, comments = comments, tabix_index = FALSE)

logging::loginfo(str_interp("Analysis completed, with results located at ${script_args$output_dir}"))




