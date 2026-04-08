#' Plot Learning / Memory Index (MI) from analearn()
#'
#' @param alr Result returned by analearn()
#' @param condition Optional trained condition to plot
#' @param ymax Upper y-axis limit (default 0.5)
#' @param ylab Y-axis label (default "MI")
#' @param label_style How to display p-values: "stars", "numeric", or "both"
#'
#' @return ggplot object
#' @export
ccPlotMI <- function(
    alr,
    condition = NULL,
    ymax = 0.5,
    ylab = "MI",
    label_style = c("stars", "numeric", "both")
) {

  label_style <- match.arg(label_style)

  if (is.null(alr$data)) {
    stop("ccPlotMI() expects output from analearn() containing $data")
  }

  df <- alr$data

  if (!is.numeric(df$p)) {
    stop(
      "ccPlotMI(): alr$data$p must be numeric.\n",
      "Re-run analearn() after changing p-value formatting."
    )
  }

  # Select condition if needed
  if (!is.null(condition)) {
    df <- df[df$condition == condition, , drop = FALSE]
  } else if ("condition" %in% names(df) && length(unique(df$condition)) > 1) {
    df <- df[df$condition == unique(df$condition)[1], , drop = FALSE]
  }

  # Preserve genotype order
  df$genotype <- factor(df$genotype, levels = df$genotype)
  control <- levels(df$genotype)[1]

  # Format labels
  p_to_stars <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p < 1e-4) "****"
    else if (p < 1e-3) "***"
    else if (p < 1e-2) "**"
    else if (p < 0.05) "*"
    else "ns"
  }

  df$label <- switch(
    label_style,
    stars   = vapply(df$p, p_to_stars, character(1)),
    numeric = format.pval(df$p, digits = 2),
    both    = paste0(format.pval(df$p, digits = 2), " ", vapply(df$p, p_to_stars, character(1)))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = genotype, y = LI)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = ylab) +
    ggplot2::coord_cartesian(ylim = c(NA, ymax)) +
    ggplot2::geom_text(
      data = df[df$genotype != control, , drop = FALSE],
      ggplot2::aes(label = label),
      nudge_y = 0.02,
      na.rm = TRUE
    )
}
