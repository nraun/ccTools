#' @keywords internal
give.n <- function(x){
  c(y = -0.05, label = length(x))
}

#' Courtship boxplots
#'
#' Plot courtship conditioning results as boxplots with pairwise
#' Wilcoxon Rank Sum tests and user-controlled p-value labels.
#'
#' @param dat A data frame generated from loadData() followed by findOutlier()
#' @param label_style How to display p-values: "stars", "numeric", or "both"
#'
#' @return A ggplot object
#' @export
ccBoxplots <- function(
    dat,
    label_style = c("stars", "numeric", "both")
) {

  label_style <- match.arg(label_style)

  stat <- dat %>%
    dplyr::group_by(genotype) %>%
    rstatix::wilcox_test(CI ~ condition, paired = FALSE) %>%
    rstatix::add_xy_position(x = "condition")

  p_to_stars <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p < 1e-4) "****"
    else if (p < 1e-3) "***"
    else if (p < 1e-2) "**"
    else if (p < 0.05) "*"
    else "ns"
  }

  stat$label <- switch(
    label_style,
    stars   = vapply(stat$p, p_to_stars, character(1)),
    numeric = format.pval(stat$p, digits = 2),
    both    = paste0(format.pval(stat$p, digits = 2),
                     " ",
                     vapply(stat$p, p_to_stars, character(1)))
  )

  ggplot2::ggplot(dat, ggplot2::aes(y = CI, x = condition)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_grid(cols = ggplot2::vars(genotype)) +
    ggpubr::stat_pvalue_manual(
      stat,
      label = "label",
      tip.length = 0
    ) +
    ggplot2::stat_summary(fun.data = give.n, geom = "text") +
    ggplot2::scale_y_continuous(expand = c(0, 0.1)) +
    ggplot2::labs(x = NULL)
}
