give.n <- function(x){
  return(c(y = -0.05, label = length(x))) 
}

ccBoxplots <- function(dat){
stat <- dat %>% group_by(genotype) %>% 
  rstatix::wilcox_test(CI ~ condition, paired=F) %>%
  rstatix::add_xy_position(x="condition")

ggplot2::ggplot(dat,aes(y=CI,x=condition)) + ggplot2::geom_boxplot() + ggplot2::facet_grid(cols=vars(genotype)) +
  ggpubr::stat_pvalue_manual(stat, tip.length = 0) + ggplot2::stat_summary(fun.data = give.n,geom="text") +
  ggplot2::scale_y_continuous(expand = c(0,0.1)) + ggplot2::labs(x=NULL)
}
