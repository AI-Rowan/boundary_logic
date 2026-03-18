############################################################
# Phase 3: Shapley values and sparse counterfactuals
# Adapted from: scripts/1.2 Run_Shapley_optimized.R
#               scripts/4.4 Local Interpretation Shapley.R
############################################################


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

#' Build a data frame of mixed start/end rows from binary masks
#'
#' Each row of \code{masks} is a binary vector: 0 means take the value from
#' \code{start}, 1 means take the value from \code{end}.
#'
#' @param start     Numeric vector length p (observed values).
#' @param end       Numeric vector length p (counterfactual values).
#' @param masks     Integer matrix [R x p] of 0/1 values.
#' @param var_names Character vector of column names length p.
#' @return Data frame with \code{var_names} columns and \code{nrow(masks)} rows.
#' @noRd
.build_rows_from_masks <- function(start, end, masks, var_names) {
  X <- sweep((1 - masks), 2, start, `*`) + sweep(masks, 2, end, `*`)
  X <- as.data.frame(X)
  colnames(X) <- var_names
  X
}


#' Exact Shapley values for a single (start, end) pair
#'
#' Uses the exact weighted marginal contribution formula for all \eqn{2^p}
#' subsets. Suitable for \eqn{p \le 14} variables; beyond that use
#' \code{.shapley_perm_one}.
#'
#' @param start      Numeric vector length p (observed).
#' @param end        Numeric vector length p (counterfactual).
#' @param model      Fitted model object.
#' @param model_type Character string identifying the model type.
#' @param rounding   Integer; rounding precision for predictions.
#' @param var_names  Character vector of variable names length p.
#' @return Named numeric vector of Shapley values length p.
#' @noRd
.shapley_exact_one <- function(start, end, model, model_type, rounding, var_names) {
  p    <- length(start)
  fact <- factorial(0:p)
  wfun <- function(k) (fact[k + 1L] * fact[p - k]) / fact[p + 1L]
  shap <- numeric(p)

  for (k in 0:(p - 1L)) {
    comb <- if (k == 0L) matrix(integer(0), nrow = 0L, ncol = 1L) else combn(p, k)
    masks_S <- if (k == 0L) {
      matrix(0L, nrow = 1L, ncol = p)
    } else {
      m_s <- matrix(0L, nrow = ncol(comb), ncol = p)
      for (cc in seq_len(ncol(comb))) m_s[cc, comb[, cc]] <- 1L
      m_s
    }
    X_S <- .build_rows_from_masks(start, end, masks_S, var_names)
    f_S <- .pred_function(model_use = model, model_type = model_type,
                          rounding = rounding, new_data = X_S)
    pairs <- list(); idx_S <- list(); idx_i <- list(); cnt <- 1L
    for (rowS in seq_len(nrow(masks_S))) {
      mask_row <- masks_S[rowS, ]
      zeros    <- which(mask_row == 0L)
      for (ii in zeros) {
        new_mask       <- mask_row; new_mask[ii] <- 1L
        pairs[[cnt]]   <- new_mask
        idx_S[[cnt]]   <- rowS
        idx_i[[cnt]]   <- ii
        cnt            <- cnt + 1L
      }
    }
    if (length(pairs) == 0L) next
    masks_Si <- do.call(rbind, pairs)
    X_Si     <- .build_rows_from_masks(start, end, masks_Si, var_names)
    f_Si     <- .pred_function(model_use = model, model_type = model_type,
                               rounding = rounding, new_data = X_Si)
    wk <- wfun(k)
    for (mm in seq_len(nrow(masks_Si))) {
      inc               <- f_Si[mm] - f_S[idx_S[[mm]]]
      shap[idx_i[[mm]]] <- shap[idx_i[[mm]]] + wk * inc
    }
  }
  shap
}


#' Approximate Shapley values via permutation sampling
#'
#' Estimates Shapley values using random permutations of the variable order.
#' Use when \eqn{p > 14}.
#'
#' @param start      Numeric vector length p (observed).
#' @param end        Numeric vector length p (counterfactual).
#' @param model      Fitted model object.
#' @param model_type Character string identifying the model type.
#' @param rounding   Integer; rounding precision for predictions.
#' @param var_names  Character vector of variable names length p.
#' @param M          Integer; number of permutations. Default \code{2048L}.
#' @param seed       Random seed. Default \code{1L}.
#' @return Named numeric vector of approximate Shapley values length p.
#' @noRd
.shapley_perm_one <- function(start, end, model, model_type, rounding, var_names,
                              M = 2048L, seed = 1L) {
  set.seed(seed)
  p        <- length(start)
  shap_sum <- numeric(p)
  for (mm in seq_len(M)) {
    ord   <- sample.int(p, p, replace = FALSE)
    masks <- matrix(0L, nrow = p + 1L, ncol = p)
    for (ii in seq_len(p)) masks[(ii + 1L), ord[ii]] <- 1L
    masks  <- apply(masks, 2L, cumsum)
    X_path <- .build_rows_from_masks(start, end, masks, var_names)
    vals   <- .pred_function(model_use = model, model_type = model_type,
                             rounding = rounding, new_data = X_path)
    d             <- diff(vals)
    shap_sum[ord] <- shap_sum[ord] + d
  }
  shap_sum / M
}


# ---------------------------------------------------------------------------
# Exported functions
# ---------------------------------------------------------------------------

#' Compute Shapley variable contributions for a local counterfactual
#'
#' Explains \emph{which variables drive the prediction from the observed point
#' to the boundary counterfactual} using Shapley values. For \eqn{p \le}
#' \code{exact_max_vars} variables, the exact formula is used; otherwise a
#' permutation approximation is applied.
#'
#' The \code{Contribute} column in the returned data frame classifies each
#' variable as \code{"Supports"} (pushes the prediction toward the boundary)
#' or \code{"Contradicts"} (pushes it away from the boundary).
#'
#' @param bl_local_result A \code{"bl_local_result"} object produced by
#'   \code{\link{bl_find_local_cf}}.
#' @param exact_max_vars  Integer; use exact computation when the number of
#'   variables is at most this value. Default \code{14L}.
#' @param approx_perm     Integer; number of permutations for the approximate
#'   method. Default \code{2048L}.
#' @param seed            Random seed for the approximate method. Default
#'   \code{1L}.
#'
#' @return A list of class \code{"bl_shapley"} with the following fields:
#'   \describe{
#'     \item{shapley_df}{Data frame; one row per variable, sorted by
#'       contribution, with columns \code{varnames}, \code{pred_data},
#'       \code{data_to_boundary}, \code{shapley_cause}, \code{pred_use},
#'       \code{Contribute}, \code{varnames_p}.}
#'     \item{pred_prob}{Scalar predicted probability of the target.}
#'     \item{pred_class}{Integer predicted class of the target.}
#'     \item{pred_boundary}{Scalar predicted probability at the boundary.}
#'     \item{row_id}{Integer row ID, or \code{NA_integer_} for external data.}
#'     \item{bl_local_result}{The \code{"bl_local_result"} passed in.}
#'   }
#'
#' @importFrom dplyr case_when
#' @export
bl_shapley <- function(bl_local_result, exact_max_vars = 14L,
                       approx_perm = 2048L, seed = 1L) {
  if (!inherits(bl_local_result, "bl_local_result"))
    stop("'bl_local_result' must be a 'bl_local_result' object.", call. = FALSE)
  if (!bl_local_result$solution_found)
    stop("No solution found in 'bl_local_result'. Run bl_find_local_cf() first.",
         call. = FALSE)

  bl_result  <- bl_local_result$bl_result
  bl_target  <- bl_local_result$bl_target
  var_names  <- bl_result$var_names
  p          <- length(var_names)

  start <- as.numeric(bl_target$x_obs[, var_names])
  end   <- as.numeric(bl_local_result$B_x[, var_names])

  if (p <= exact_max_vars) {
    shap <- .shapley_exact_one(start, end, bl_result$model, bl_result$model_type,
                               bl_result$rounding, var_names)
  } else {
    shap <- .shapley_perm_one(start, end, bl_result$model, bl_result$model_type,
                              bl_result$rounding, var_names,
                              M = as.integer(approx_perm), seed = seed)
  }
  names(shap) <- var_names

  pred_prob     <- bl_target$pred_prob
  pred_class    <- bl_target$pred_class
  pred_use      <- pred_class
  pred_boundary <- bl_local_result$B_pred
  data_to_boundary <- end - start

  # Ordering: class 0 → decreasing (most negative first); class 1 → ascending
  dec_order  <- (pred_use == 1L)
  dist_order <- order(shap, decreasing = dec_order)

  df <- data.frame(
    varnames         = var_names,
    pred_data        = start,
    data_to_boundary = data_to_boundary,
    shapley_cause    = shap,
    pred_use         = pred_use,
    stringsAsFactors = FALSE
  )
  df$Contribute <- dplyr::case_when(
    df$pred_use == 0L & df$shapley_cause >= 0 ~ "Supports",
    df$pred_use == 0L & df$shapley_cause <  0 ~ "Contradicts",
    df$pred_use == 1L & df$shapley_cause >= 0 ~ "Contradicts",
    df$pred_use == 1L & df$shapley_cause <  0 ~ "Supports",
    TRUE ~ "Unknown"
  )

  df <- df[dist_order, , drop = FALSE]
  df$varnames   <- factor(df$varnames, levels = df$varnames)
  df$varnames_p <- paste0(
    df$varnames, ": ",
    round(df$pred_data, 2), " -> ", round(df$data_to_boundary, 2)
  )
  df$varnames_p <- factor(df$varnames_p, levels = df$varnames_p)

  structure(
    list(
      shapley_df      = df,
      pred_prob       = pred_prob,
      pred_class      = pred_class,
      pred_boundary   = pred_boundary,
      row_id          = bl_target$row_id,
      bl_local_result = bl_local_result
    ),
    class = "bl_shapley"
  )
}


#' Plot Shapley contributions for a local counterfactual
#'
#' Returns a \pkg{ggplot2} bar chart showing each variable's Shapley
#' contribution to the path from the observed point to the boundary.
#' Variables are coloured by whether they \emph{support} or
#' \emph{contradict} the prediction direction.
#'
#' @param x   A \code{"bl_shapley"} object from \code{\link{bl_shapley}}.
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot} object (invisible printing — assign or call
#'   explicitly to display).
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_label scale_fill_manual
#'   scale_color_manual theme_bw theme labs
#' @export
plot.bl_shapley <- function(x, ...) {
  df        <- x$shapley_df
  row_label <- if (is.na(x$row_id)) "external" else as.character(x$row_id)

  ggplot2::ggplot(
    df,
    ggplot2::aes(fill = Contribute, colour = Contribute,
                 x = shapley_cause, y = varnames_p)
  ) +
    ggplot2::scale_fill_manual(
      name   = "Contribution",
      values = c("Supports" = "darkblue", "Contradicts" = "grey70",
                 "Unknown"  = "lightgrey")
    ) +
    ggplot2::scale_color_manual(
      name   = "Contribution",
      values = c("Supports" = "darkblue", "Contradicts" = "grey70",
                 "Unknown"  = "lightgrey")
    ) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::geom_label(
      ggplot2::aes(label = round(shapley_cause, 3L), x = shapley_cause),
      colour = "white", show.legend = FALSE
    ) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x        = "[-> change to reach boundary]: Impact on prediction",
      y        = "Variable (observed value -> change required)",
      title    = "Shapley Contribution Plot",
      subtitle = sprintf(
        "ID: %s | Class: %d | Pred: %.3f | Boundary pred: %.3f",
        row_label, x$pred_class, x$pred_prob, x$pred_boundary
      )
    )
}


#' Print method for bl_shapley
#'
#' @param x A \code{"bl_shapley"} object.
#' @param ... Currently ignored.
#' @export
print.bl_shapley <- function(x, ...) {
  cat("-- bl_shapley --\n")
  rid <- if (is.na(x$row_id)) "external" else as.character(x$row_id)
  cat(sprintf("  Row ID         : %s\n",  rid))
  cat(sprintf("  Pred prob      : %.4f\n", x$pred_prob))
  cat(sprintf("  Pred class     : %d\n",   x$pred_class))
  cat(sprintf("  Boundary pred  : %.4f\n", x$pred_boundary))
  cat("\n  Shapley table:\n")
  df_print <- x$shapley_df[, c("varnames", "pred_data", "data_to_boundary",
                                "shapley_cause", "Contribute"),
                            drop = FALSE]
  df_print$varnames         <- as.character(df_print$varnames)
  df_print$pred_data        <- round(df_print$pred_data, 4L)
  df_print$data_to_boundary <- round(df_print$data_to_boundary, 4L)
  df_print$shapley_cause    <- round(df_print$shapley_cause, 4L)
  print(df_print, row.names = FALSE)
  invisible(x)
}


#' Create a sparse counterfactual by zeroing contradicting variables
#'
#' Builds a sparse version of the boundary counterfactual by retaining only
#' the variables whose Shapley contribution \emph{supports} the prediction
#' direction. Variables classified as \code{"Contradicts"} or \code{"Unknown"}
#' revert to their observed values. Additionally, any variable constrained as
#' \code{"fixed"} in \code{\link{set_filters}} always reverts to its observed
#' value regardless of its Shapley classification, because the \code{"fixed"}
#' constraint permits only a \eqn{\pm 0.5} search window and the observed value
#' is the intended sparse target.
#'
#' @param bl_shapley_result A \code{"bl_shapley"} object from
#'   \code{\link{bl_shapley}}.
#' @param round_to Positive numeric or \code{NULL}. If not \code{NULL}, rounds
#'   supporting variable values to the nearest multiple of this value (e.g.
#'   \code{round_to = 0.5} rounds to the nearest 0.5; \code{round_to = 1}
#'   rounds to the nearest integer). Default \code{NULL} (no rounding).
#'
#' @return A list of class \code{"bl_sparse_result"} with the following fields:
#'   \describe{
#'     \item{x_sparse}{Single-row data frame; the sparse counterfactual in
#'       X-space.}
#'     \item{pred_sparse}{Scalar predicted probability of the sparse CF.}
#'     \item{solution_valid}{Logical; \code{TRUE} if the sparse CF flips the
#'       predicted class.}
#'     \item{shapley_df}{The annotated Shapley data frame with columns
#'       \code{used_in_sparse} and \code{x_sparse_val}.}
#'     \item{bl_shapley}{The \code{"bl_shapley"} object passed in.}
#'   }
#'
#' @export
bl_find_sparse_cf <- function(bl_shapley_result, round_to = NULL) {
  if (!inherits(bl_shapley_result, "bl_shapley"))
    stop("'bl_shapley_result' must be a 'bl_shapley' object from bl_shapley().",
         call. = FALSE)

  df        <- bl_shapley_result$shapley_df
  bl_local  <- bl_shapley_result$bl_local_result
  bl_result <- bl_local$bl_result
  bl_target <- bl_local$bl_target
  var_names <- bl_result$var_names

  x_obs <- as.numeric(bl_target$x_obs[, var_names])
  B_x   <- as.numeric(bl_local$B_x[, var_names])
  names(x_obs) <- var_names
  names(B_x)   <- var_names

  # Build sparse CF: use CF value only for "Supports" variables;
  # "Contradicts" and "Unknown" revert to observed value.
  x_sparse <- x_obs
  for (nm in var_names) {
    row_nm <- df[as.character(df$varnames) == nm, , drop = FALSE]
    if (nrow(row_nm) == 0L) next
    if (row_nm$Contribute[1L] == "Supports") {
      new_val <- B_x[[nm]]
      if (!is.null(round_to)) new_val <- round(new_val / round_to) * round_to
      x_sparse[[nm]] <- new_val
    }
  }

  # "fixed"-constrained variables always revert to observed value in sparse CF,
  # regardless of Shapley classification (set_filters stored in bl_local_result).
  filters <- bl_local$set_filters
  if (!is.null(filters) && length(filters$constraints) > 0L) {
    for (nm in names(filters$constraints)) {
      v <- filters$constraints[[nm]]
      if (is.character(v) && v == "fixed") {
        x_sparse[[nm]] <- x_obs[[nm]]
      }
    }
  }

  x_sparse_df <- as.data.frame(t(x_sparse))
  colnames(x_sparse_df) <- var_names

  pred_sparse <- .pred_function(
    model_use  = bl_result$model,
    model_type = bl_result$model_type,
    rounding   = bl_result$rounding,
    new_data   = x_sparse_df
  )

  pred_class_target <- bl_target$pred_class
  pred_class_sparse <- as.integer(pred_sparse[[1L]] >= bl_result$cutoff)
  solution_valid    <- (pred_class_sparse != pred_class_target)

  # Annotate data frame
  df_out                <- df
  df_out$used_in_sparse <- df_out$Contribute == "Supports"
  df_out$x_sparse_val   <- x_sparse[as.character(df_out$varnames)]

  structure(
    list(
      x_sparse       = x_sparse_df,
      pred_sparse    = pred_sparse[[1L]],
      solution_valid = solution_valid,
      shapley_df     = df_out,
      bl_shapley     = bl_shapley_result
    ),
    class = "bl_sparse_result"
  )
}


#' Plot a sparse counterfactual result
#'
#' Renders the rotated local biplot via \code{\link{plot.bl_local_result}}
#' (which uses the biplotEZ pipeline) and then overlays the sparse
#' counterfactual as a second cross marker. A blue cross indicates the sparse
#' CF successfully flips the predicted class; red indicates it does not.
#' All \code{\link{plot.bl_local_result}} parameters are forwarded through
#' \code{...}; the local summary is suppressed and replaced by a combined
#' sparse CF summary.
#'
#' @param x         A \code{"bl_sparse_result"} object from
#'   \code{\link{bl_find_sparse_cf}}.
#' @param show_arrows Logical; passed to \code{\link{plot.bl_local_result}}
#'   to control the CF arrow. Default \code{TRUE}.
#' @param arrow_col Character; colour for the full CF cross and arrow. Default
#'   \code{"grey30"}.
#' @param ...       Additional arguments passed to
#'   \code{\link{plot.bl_local_result}} (e.g. \code{no_grid},
#'   \code{label_offset_var}, \code{which}).
#'
#' @importFrom graphics points legend
#' @export
plot.bl_sparse_result <- function(x,
                                  show_arrows = TRUE,
                                  arrow_col   = "grey30",
                                  ...) {
  # Render underlying local biplot; suppress its summary (we print our own)
  plot(x$bl_shapley$bl_local_result,
       show_arrows   = show_arrows,
       arrow_col     = arrow_col,
       print_summary = FALSE,
       ...)

  bl_local    <- x$bl_shapley$bl_local_result
  bl_result   <- bl_local$bl_result
  var_names   <- bl_result$var_names
  X_center    <- bl_result$X_center
  X_sd        <- bl_result$X_sd
  standardise <- bl_result$standardise
  Vr_rot      <- bl_local$Vr_rot

  # Project sparse CF into rotated Z-space
  x_sparse_vec <- as.numeric(x$x_sparse[, var_names])
  sv           <- if (isTRUE(standardise)) X_sd else rep(1, length(var_names))
  x_sparse_st  <- (x_sparse_vec - X_center) / sv
  z_sparse     <- matrix(x_sparse_st %*% Vr_rot, nrow = 1L)

  pt_col <- if (x$solution_valid) "green3" else "yellow2"

  graphics::points(z_sparse[, 1L], z_sparse[, 2L],
                   pch = 4L, col = pt_col, cex = 1.4, lwd = 2)

  graphics::legend(
    "topright",
    legend = c("Full CF",
               if (x$solution_valid) "Sparse CF (valid)" else "Sparse CF (invalid)"),
    col    = c(arrow_col, pt_col),
    pch    = 4L,
    pt.lwd = c(1.5, 2.0),
    pt.cex = c(1.2, 1.4),
    bty    = "n",
    cex    = 0.8
  )

  # Sparse CF summary
  bl_target <- bl_local$bl_target
  row_label <- if (is.na(bl_target$row_id)) "external" else
                 as.character(bl_target$row_id)
  cat("\n--- Sparse CF summary ---\n")
  cat(sprintf("  Target         : %s\n",         row_label))
  cat(sprintf("  Observed pred  : %.4f  (class %d)\n",
              x$bl_shapley$pred_prob, x$bl_shapley$pred_class))
  cat(sprintf("  Full CF pred   : %.4f\n",       x$bl_shapley$pred_boundary))
  cat(sprintf("  Sparse CF pred : %.4f\n",       x$pred_sparse))
  cat(sprintf("  Solution valid : %s\n",         x$solution_valid))
  cat("-------------------------\n")

  invisible(x)
}


#' Print method for bl_sparse_result
#'
#' @param x A \code{"bl_sparse_result"} object.
#' @param ... Currently ignored.
#' @export
print.bl_sparse_result <- function(x, ...) {
  cat("-- bl_sparse_result --\n")
  cat(sprintf("  Solution valid : %s\n", x$solution_valid))
  cat("\n  Predictions:\n")
  cat(sprintf("    Observed     : %.4f  (class %d)\n",
              x$bl_shapley$pred_prob, x$bl_shapley$pred_class))
  cat(sprintf("    Full CF      : %.4f\n", x$bl_shapley$pred_boundary))
  cat(sprintf("    Sparse CF    : %.4f\n", x$pred_sparse))
  cat("\n  Variable summary:\n")
  df <- x$shapley_df
  bl_local  <- x$bl_shapley$bl_local_result
  bl_target <- bl_local$bl_target
  var_names <- bl_local$bl_result$var_names

  x_obs_vec <- as.numeric(bl_target$x_obs[, var_names])
  B_x_vec   <- as.numeric(bl_local$B_x[, var_names])

  df_print <- data.frame(
    variable      = as.character(df$varnames),
    x_obs         = round(x_obs_vec[match(as.character(df$varnames), var_names)], 4L),
    B_x           = round(B_x_vec[match(as.character(df$varnames), var_names)], 4L),
    x_sparse      = round(df$x_sparse_val, 4L),
    used_in_sparse = df$used_in_sparse,
    stringsAsFactors = FALSE
  )
  print(df_print, row.names = FALSE)
  invisible(x)
}
