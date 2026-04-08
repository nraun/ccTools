#' Analyze learning index (LI) and genotype differences via bootstrap
#'
#' Computes learning index (LI) per genotype and trained condition:
#' LI = (CInaive - CItrained) / CInaive
#'
#' If controlGenotype is provided, computes bootstrap confidence intervals
#' and two-sided p-values for the difference in LI versus control:
#' (LI_control - LI_genotype).
#'
#' @param data Data frame with columns: genotype, condition, CI.
#' @param nboot Number of bootstrap resamples (default 10000).
#' @param naivelevel Which condition level is "naive" baseline (default "N").
#' @param seed Optional seed for reproducibility (default NA).
#' @param controlGenotype Reference genotype for inference (default NULL).
#' @param method Aggregation function for CI within genotype/condition: "mean" or "median".
#' @param conf Confidence level for bootstrap CI (default 0.95).
#'
#' @return A list with vectors genotype, LI, p and full results data.frame in $data.
#' @export
analearn <- function(
    data,
    nboot = 10000,
    naivelevel = "N",
    seed = NA,
    controlGenotype = NULL,
    method = c("mean", "median"),
    conf = 0.95
) {
  if (is.null(controlGenotype)) {
    controlGenotype <- attr(data, "refgenotype")
  }

  if (is.null(controlGenotype)) {
    stop(
      "controlGenotype not provided and not found in data attributes.\n",
      "Either:\n",
      " - call loadData(..., refgenotype = 'X') first, or\n",
      " - pass controlGenotype explicitly to analearn()."
    )
  }

  method <- match.arg(method)

  req <- c("genotype", "condition", "CI")
  miss <- setdiff(req, names(data))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

  data <- as.data.frame(data)
  if (!is.na(seed)) set.seed(seed)

  # Factors + relevel
  data$condition <- as.factor(data$condition)
  if (!naivelevel %in% levels(data$condition)) {
    stop("naivelevel='", naivelevel, "' not found in condition levels: ",
         paste(levels(data$condition), collapse = ", "))
  }
  data$condition <- stats::relevel(data$condition, naivelevel)

  data$genotype <- as.factor(data$genotype)
  if (!is.null(controlGenotype)) {
    if (!controlGenotype %in% levels(data$genotype)) {
      stop("controlGenotype='", controlGenotype, "' not found in genotype levels: ",
           paste(levels(data$genotype), collapse = ", "))
    }
    data$genotype <- stats::relevel(data$genotype, ref = controlGenotype)
  }

  stat_fun <- if (method == "mean") base::mean else stats::median

  # Aggregate CI by genotype/condition
  aggr <- stats::aggregate(CI ~ genotype + condition, data = data, FUN = stat_fun)

  # Build LI table
  naive_ci <- aggr[aggr$condition == naivelevel, c("genotype", "CI")]
  names(naive_ci)[2] <- "CInaive"

  trained <- aggr[aggr$condition != naivelevel, c("genotype", "condition", "CI")]
  names(trained)[3] <- "CItrained"

  res <- merge(trained, naive_ci, by = "genotype", all.x = TRUE)
  res$LI <- (res$CInaive - res$CItrained) / res$CInaive

  # Initialize inference columns
  res$LIdif  <- NA_real_
  res$LL95CI <- NA_real_
  res$UL95CI <- NA_real_
  res$p <- NA_real_

  # Bootstrap diffs vs control (optional)
  if (!is.null(controlGenotype) && nboot > 0) {
    control <- levels(data$genotype)[1]
    geno_levels <- levels(data$genotype)
    cond_levels <- levels(data$condition)
    trained_levels <- setdiff(cond_levels, naivelevel)
    mutants <- setdiff(geno_levels, control)

    # Split row indices by genotype/condition
    grp <- split(seq_len(nrow(data)), interaction(data$genotype, data$condition, drop = TRUE))
    grp_sizes <- vapply(grp, length, integer(1))
    grp_names <- do.call(rbind, strsplit(names(grp), "\\."))
    grp_keys <- data.frame(genotype = grp_names[, 1], condition = grp_names[, 2], stringsAsFactors = FALSE)

    # Bootstrapped diffs: (LI_control - LI_mutant) per trained condition
    ncol_boot <- length(trained_levels) * length(mutants)
    boot <- matrix(NA_real_, nrow = nboot, ncol = ncol_boot)
    colnames(boot) <- as.vector(outer(trained_levels, mutants, paste, sep = "__"))

    for (b in seq_len(nboot)) {
      ci_star <- mapply(function(idx, n) stat_fun(data$CI[sample(idx, n, replace = TRUE)]),
                        idx = grp, n = grp_sizes, SIMPLIFY = TRUE)

      ag <- cbind(grp_keys, CI = as.numeric(ci_star))
      ag$genotype <- factor(ag$genotype, levels = geno_levels)
      ag$condition <- factor(ag$condition, levels = cond_levels)

      ag_naive <- ag[ag$condition == naivelevel, c("genotype", "CI")]
      names(ag_naive)[2] <- "CInaive"

      ag_tr <- ag[ag$condition != naivelevel, c("genotype", "condition", "CI")]
      names(ag_tr)[3] <- "CItrained"

      tmp <- merge(ag_tr, ag_naive, by = "genotype", all.x = TRUE)
      tmp$LI <- (tmp$CInaive - tmp$CItrained) / tmp$CInaive

      k <- 0
      for (cond in trained_levels) {
        li_c <- tmp$LI[tmp$genotype == control & tmp$condition == cond]
        for (mut in mutants) {
          k <- k + 1
          li_m <- tmp$LI[tmp$genotype == mut & tmp$condition == cond]
          boot[b, k] <- li_c - li_m
        }
      }
    }

    alpha <- (1 - conf) / 2
    eps <- 1 / nboot

    for (cond in trained_levels) {
      li_control_obs <- res$LI[res$genotype == control & res$condition == cond]

      for (mut in mutants) {
        key <- paste(cond, mut, sep = "__")
        dist <- boot[, key]

        q <- stats::quantile(dist, probs = c(alpha, 1 - alpha), na.rm = TRUE)
        p1 <- mean(dist > 0, na.rm = TRUE)
        p  <- 2 * min(p1, 1 - p1)
        if (p == 0) p <- eps

        irow <- which(res$genotype == mut & res$condition == cond)
        if (length(irow) == 1) {
          res$LIdif[irow]  <- li_control_obs - res$LI[irow]
          res$LL95CI[irow] <- q[[1]]
          res$UL95CI[irow] <- q[[2]]

          res$p[irow] <- p
        }
      }
    }
  }

  # Stable ordering
  res <- res[order(res$condition, res$genotype), ]

  # Backward-friendly: if only one trained condition, return vectors over genotype
  trained_levels <- setdiff(levels(data$condition), naivelevel)
  if (length(trained_levels) == 1) {
    out <- res[res$condition == trained_levels[1], ]
    return(structure(list(
      genotype = as.character(out$genotype),
      LI = out$LI,
      p = out$p,
      data = res
    ), class = "analearn_result"))
  }

  structure(list(
    genotype = as.character(res$genotype),
    learningcondition = as.character(res$condition),
    LI = res$LI,
    p = res$p,
    data = res
  ), class = "analearn_result")
}
