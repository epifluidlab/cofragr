#!/usr/bin/env Rscript

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

library(here)

# example
script_args <- list(
  input = here::here("cofrag/results/20211108-week-45/cofrag/EE86217.GRCh38.cofrag_cm.bed.gz"),
  output_dir = here::here("sandbox/"),
  sample_id = "IH01",
  res = 500e3L,
  standard_compartment = here("cofrag/results/20211004-week-40/wbc/compartment/wbc_rep1.hg38.comps.juicer.bedGraph"),
  chroms = c("20", "21", "22"),
  method = c("juicer", "lieberman"),
  smooth = c(1L, 3L),
  genome = "GRCh38",
  java = "java"
)


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
    assert_that(is_scalar_character(genome) && genome %in% c("GRCh38", "GRCh37", "hs37-1kg"))
    assert_that(length(method) >= 1 &&
                  all(method %in% c("juicer", "lieberman", "obs_exp")))
    # assert_that(length(smooth) >= 1 && is_integer(smooth) && all(smooth > 0))
  }) %>% invisible()
} 


if (interactive()) {
  library(tidyverse)
  library(magrittr)
  library(rlang)
  library(assertthat)

  if (is.null(get0("script_args")))
    stop_quietly()
} else {
  library(rlang)
  library(assertthat)
  
  # Run in CLI script mode
  parser <- optparse::OptionParser(
    option_list = list(
      optparse::make_option(c("-i", "--input"), type = "character", help = "Input files (BEDPE format)"),
      optparse::make_option(c("-o", "--output-dir", type = "character", help = "Directory for output")),
      optparse::make_option(c("-s", "--sample-id"), type = "character", help = "Sample name. This affects output file names"),
      optparse::make_option(c("--res"), type = "integer", default = 500e3L, help = "Resolution used in the compartment analysis [500000]"),
      optparse::make_option(c("--chroms"), default = NULL,
                            help = "Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used [NULL]"),
      optparse::make_option(c("--method"), default = "juicer:lieberman", help = "Methods for compartment analysis. If multiple methods are involved, must be separated by colons [juicer:lieberman]"),
      optparse::make_option(c("--smooth"), default = "1:3", help = "Apply moving average on the compartment track. Sometimes this can remove some quirks and make the data more similar to Hi-C compartment tracks [1:3]"),
      optparse::make_option(
        c("--standard-compartment"),
        default = NULL,
        help = "A standard compartment track (BED format) for calculating compartment correlation scores. Can be `gc`, `gene.density` or a file path [NULL]"
      ),
      optparse::make_option(
        c("-g", "--genome"),
        default = NULL,
        help = "Reference genome of the dataset [GRCh38]",
      ),
      # optparse::make_option(c("--juicer"), default = NULL,
      #                       help = "Path to the .jar file Juicer tools. If not provided, the one shipped with hictools package will be used"),
      optparse::make_option(c("--java"), default = "java",
                            help = "Path to java [java]")
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

if (script_args$standard_compartment == "gc") {
  standard_comp <- local({
    env <- new.env()
    res <- paste0(script_args$res %/% 1000, "kbp")
    track_name <- stringr::str_interp("gc.GRCh37.${res}")
    data(list = track_name,
         package = "hictools",
         envir = env)
    env[[track_name]]
  })
} else if (!is_null(script_args$standard_compartment)) {
  standard_comp <-
    bedtorch::read_bed(file_path = script_args$standard_compartment,
                       genome = script_args$genome)
} else
  standard_comp <- NULL


cofrag_comp_results <- script_args$method %>%
  # pmap(function(chrom, method, smooth, metric) {
  map(function(method) {
    metric <- "score"
    
    # logging::loginfo(str_interp("Processing: chrom:${chrom} ${method} smooth:${smooth} ${metric}"))
    logging::loginfo(str_interp("Processing using ${method} ..."))
    
    hic_data <- bedtorch::read_bed(script_args$input, genome = script_args$genome)
    
    if ("bootstrap" %in% colnames(GenomicRanges::mcols(hic_data))) {
      hic_file <- tempfile(fileext = ".bed")
      on.exit(file.remove(hic_file), add = TRUE)
      
      hic_data <- hic_data[hic_data$bootstrap == 1]
      bedtorch::write_bed(hic_data, file_path = hic_file)
    } else {
      hic_file <- script_args$input
    }
    
    hic_data <-
      hictools::load_hic_genbed(
        hic_file,
        # chrom,
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
                                chrom = all_chroms,
                                standard = standard_comp,
                                oe = "juicer",
                                genome = script_args$genome,
                                java = script_args$java)
    comps$smooth <- 1L
    
    comps_smoothed <- setdiff(script_args$smooth, 1) %>%
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
    comps$method <- method
    comps$metric <- metric
    
    comps$score <- local({
      score <- comps$score
      ifelse(is.infinite(score) | is.na(score), NA, score)
    })
    
    # if (!is_null(standard_comp))
    #   correlation <- hictools::comp_correlation(comps, standard_comp)
    # else
    #   correlation <- NULL
    
    comps
    # list(comps = comps, correlation = correlation)
  }) 

comps <- do.call(c, cofrag_comp_results)


comps[comps$method == "juicer" & comps$smooth == 1] %>%
  hictools::plot_compartment(chrom = "20")

comps %>% bedtorch::as.bedtorch_table() %>% as_tibble() %>%
  filter(smooth == 1) %>%
  mutate(pos = paste0(chrom, ":", start)) %>%
  select(pos, score, method) %>%
  pivot_wider(names_from = "method", values_from = "score") %>%
  ggpubr::ggscatter(x = "juicer", y = "lieberman")

# %>%
#   do.call(what = c, args = .)

# comps <- cofrag_comp_results %>% map(~ .$comps) %>% do.call(what = c, args = .)

# cofrag_correlation <-
#   map2_dfr(seq.int(nrow(cofrag_comp_grid)), 
#            cofrag_comp_results, 
#            function(grid_idx, result) {
#     c(cofrag_comp_grid[grid_idx,], result$correlation)
#   })

system(str_interp("mkdir -p ${script_args$output_dir}"))
cofrag_comp_file <- str_interp("${script_args$output_dir}/${script_args$sample_id}.compartment.bed.gz")
cofrag_comp_pdf <- str_interp("${script_args$output_dir}/${script_args$sample_id}.compartment.pdf")
cofrag_correlation_file <- str_interp("${script_args$output_dir}/${script_args$sample_id}.compartment_correlation.tsv")

logging::loginfo(str_interp("Write compartment scores: ${cofrag_comp_file}"))
logging::loginfo(str_interp("Write compartment correlation: ${cofrag_correlation_file}"))
bedtorch::write_bed(comps, file_path = cofrag_comp_file, comments = comments, tabix_index = FALSE)
# write_tsv(cofrag_correlation, file = cofrag_correlation_file)

logging::loginfo(str_interp("Write compartment plots: ${cofrag_comp_pdf}"))

all_chroms %>%
  map(function(chrom) {
    comps <- comps[comps$method == script_args$method[1] & comps$smooth == script_args$smooth[1]]
    
    hictools::plot_compartment(comps, chrom = chrom) +
      labs(subtitle = paste0("Chrom: ", chrom)) +
      xlab("") + ylab("") + 
      theme(legend.position = "none")
  }) %>%
  cowplot::plot_grid(plotlist = ., ncol = 1) %>%
  ggsave(
    filename = basename(cofrag_comp_pdf),
    device = "pdf",
    plot = .,
    path = dirname(cofrag_comp_pdf),
    width = 8,
    height = 2 * length(all_chroms),
    limitsize = FALSE
  )

if (!is_null(standard_comp)) {
  cofrag_vs_standard_pdf <- str_interp("${script_args$output_dir}/${script_args$sample_id}.vs_standard.pdf")
  
  logging::loginfo(str_interp("Write compartment vs. standard plots: ${cofrag_vs_standard_pdf}"))
  
  local({
    comps <-
      comps[comps$method == script_args$method[1] &
              comps$smooth == script_args$smooth[1]]
    hits <- GenomicRanges::findOverlaps(comps, standard_comp)
    x <- comps[S4Vectors::queryHits(hits)]$score
    y <-
      GenomicRanges::mcols(standard_comp[S4Vectors::subjectHits(hits)])[[1]]
    dt <- tibble(x = x, Standard = y)
    colnames(dt)[1] <- script_args$sample_id
    
    min_lim <- min(c(x, y), na.rm = TRUE)
    max_lim <- max(c(x, y), na.rm = TRUE)
    
    plot <- dt %>%
      ggpubr::ggscatter(
        x = script_args$sample_id,
        y = "Standard",
        size = 0.5,
        alpha = 0.5,
        color = "steelblue"
      ) +
      geom_abline(slope = 1,
                  linetype = "solid",
                  color = "gray") + coord_fixed() +
      xlim(min_lim, max_lim) + ylim(min_lim, max_lim) +
      ggpubr::stat_cor(label.x = min_lim, label.y = max_lim)
    
    ggsave(
      filename = basename(cofrag_vs_standard_pdf),
      device = "pdf",
      plot = plot,
      path = dirname(cofrag_vs_standard_pdf),
      width = 6,
      height = 6,
      limitsize = FALSE
    )
  })
}

logging::loginfo(str_interp("Analysis completed, with results located at ${script_args$output_dir}"))




