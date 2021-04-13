# R package for cofragmentation analysis

#' @importFrom dplyr `%>%` tibble as_tibble mutate rename select filter group_by ungroup summarize arrange inner_join lag
#' @importFrom stringr str_interp str_split str_detect str_sort str_match
#' @importFrom purrr map2_dfr map_dbl map
#' @importFrom tidyr expand_grid
#' @importFrom magrittr `%<>%`
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps mcols `mcols<-` pintersect ranges `ranges<-` width `width<-` seqnames start `start<-` end `end<-` seqinfo `seqinfo<-`
#' @importFrom S4Vectors Rle DataFrame from to
#' @importFrom parallel makeForkCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG `%dorng%`
#' @importFrom foreach foreach `%do%` `%dopar%`
NULL
