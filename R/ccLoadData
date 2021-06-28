loadData <- function(datname = NA, refgenotype = NA, naivelevel = "N"){

if (is.na(datname)) datname <- file.choose()
datfile <- read.table(datname,header = T,stringsAsFactors = FALSE)
datfile$genotype <- as.factor(datfile$genotype)
if (is.na(refgenotype) || !(refgenotype %in% levels(datfile$genotype))){
  cat('The following genotypes are observed:\n')
  for(i in seq_along(levels(datfile$genotype))) {
    cat(i,': ',levels(datfile$genotype)[i],'\n', sep = "")
  }
  cat('Please indicate the reference genotype by giving a number from 1 to ',length(levels(datfile$genotype)),' followed by <Enter>.\n',sep = "")
  refnum <- readline()
  refgenotype <- levels(datfile$genotype)[as.numeric(refnum)]}

datfile$genotype <- relevel(datfile$genotype,refgenotype)
datfile$condition <- relevel(as.factor(datfile$condition),naivelevel)
return(datfile)
}
