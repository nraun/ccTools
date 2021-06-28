gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# outliers as per https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-123
findOutlier<-function(dat, Q=0.05,outliers=F){
  con <- levels(dat$condition)
  gen <- levels(dat$genotype)
  dat2 <- NULL
  for (l in 1:length(mut)){
    for (j in 1:length(con)){
      tmp <- subset(dat,genotype ==mut[[l]] & condition ==con[[j]])
      tmp$res <- abs(lm(CI~I(seq(1:length(CI)) * 0 + gm_mean(CI)), data=tmp)$residuals)
      tmp <- tmp[order(tmp$res),]
      n <- length(tmp$res)
      rsdr <- quantile(tmp$res,0.6827)[[1]] * (n/(n-2))
      tmp$outlier <- FALSE
      for (i in round(0.7*n):n){
        a <- (Q * (n-(i-1)))/n
        t <- tmp$res[i]/rsdr
        p <- 2 * pt(t, df=n-2, lower.tail = F)
        ifelse(p<a, (tmp$outlier[i]<-TRUE),(tmp$outlier[i]<-FALSE))
      }
      dat2<-rbind(dat2,tmp)
    }
  }
  dat2$res <- NULL
  ifelse(outliers==TRUE,return(dat2),return(subset(dat2,outlier==FALSE)))
}
