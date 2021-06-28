give.n <- function(x){
  return(c(y = -0.05, label = length(x))) 
}

ccBoxplots <- function(dat){
stat <- dat %>% group_by(genotype) %>% 
  rstatix::wilcox_test(CI ~ condition, paired=F) %>%
  rstatix::add_xy_position(x="condition")

ggplot(dat,aes(y=CI,x=condition)) + geom_boxplot() + facet_grid(cols=vars(genotype)) +
  ggpubr::stat_pvalue_manual(stat, tip.length = 0) + stat_summary(fun.data = give.n,geom="text") +
  scale_y_continuous(expand = c(0,0.1)) + labs(x=NULL)
}
