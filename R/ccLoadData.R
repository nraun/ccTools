#' Load data table
#'
#' A function that loads in your analearn output file and sets the reference genotypes and conditions
#'
#' @param datname File name of tab delimited output of analearn. Leave NA to choose file in file browser.
#' @param refgenotype Character indicating the reference genotype (column name: 'genotype'). Leave NA to select interactively.
#' @param naivelevel Character indicating the naive condition (column name: 'condition'). Default is 'N'.
#' @return A data table with genotype and conditions as factors releveled to the indicated reference conditions.
#' @export

loadData <- function(datname = NA, refgenotype = NA, naivelevel = "N"){

  if (is.na(datname)) datname <- file.choose()
  datfile <- read.table(datname,header = T,stringsAsFactors = FALSE)
  datfile$genotype <- as.factor(datfile$genotype)
  if(!is.na(refgenotype) & refgenotype=="return"){
    return(levels(datfile$genotype))
  }
  if(!is.na(naivelevel) & naivelevel=="return"){
    return(levels(as.factor(datfile$condition)))
  }
  if (is.na(refgenotype) || !(refgenotype %in% levels(datfile$genotype))){
    cat('The following genotypes are observed:\n')
    for(i in seq_along(levels(datfile$genotype))) {
      cat(i,': ',levels(datfile$genotype)[i],'\n', sep = "")
    }
    cat('Please indicate the reference genotype by giving a number from 1 to ',length(levels(datfile$genotype)),' followed by <Enter>.\n',sep = "")
    refnum <- readline()
    refgenotype <- levels(datfile$genotype)[as.numeric(refnum)]
  }
  datfile$genotype <- relevel(datfile$genotype,refgenotype)
  datfile$condition <- relevel(as.factor(datfile$condition),naivelevel)
  return(datfile)
}
