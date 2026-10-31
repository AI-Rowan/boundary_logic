############################################################
# Phase 3: Shapley values and sparse counterfactuals
# Adapted from: scripts/1.2 Run_Shapley_optimized.R
#               scripts/4.4 Local Interpretation Shapley.R
############################################################

# Bare column-name symbols referenced via NSE in plot.bl_shapley()'s
# ggplot2::aes() calls -- declared so R CMD check does not flag them as
# undefined globals.
utils::globalVariables(c("Contribute", "shapley_cause", "varnames_p"))


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
#' @param masks     Integer matrix `[R x p]` of 0/1 values.
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
#' @param var_names  Character vector of variable names length p.
#' @return Named numeric vector of Shapley values length p.
#' @noRd
.shapley_exact_one <- function(start, end, model, model_type, var_names) {
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
                          new_data = X_S)
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
                               new_data = X_Si)
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
#' @param var_names  Character vector of variable names length p.
#' @param M          Integer; number of permutations. Default \code{2048L}.
#' @param seed       Random seed. \code{NULL} (default) uses the current RNG
#'   state. Supply an integer for reproducible results.
#' @return Named numeric vector of approximate Shapley values length p.
#' @noRd
.shapley_perm_one <- function(start, end, model, model_type, var_names,
                              M = 2048L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p        <- length(start)
  shap_sum <- numeric(p)
  for (mm in seq_len(M)) {
    ord   <- sample.int(p, p, replace = FALSE)
    masks <- matrix(0L, nrow = p + 1L, ncol = p)
    for (ii in seq_len(p)) masks[(ii + 1L), ord[ii]] <- 1L
    masks  <- apply(masks, 2L, cumsum)
    X_path <- .build_rows_from_masks(start, end, masks, var_names)
    vals   <- .pred_function(model_use = model, model_type = model_type,
                             new_data = X_path)
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
#' variable as \code{"Supports"} (pushes the prediction toward the boundary),
#' \code{"Contradicts"} (pushes it away from the boundary), or
#' \code{"No effect"} (Shapley value zero, up to \code{1e-8}). A zero value
#' arises when the variable's change does not move the prediction in any
#' coalition, for example a small change that crosses no split of a tree
#' model. The PhD research script counted zero as \code{"Supports"} (for a
#' class-0 target); \code{"No effect"} keeps such a variable out of the sparse
#' counterfactual (\code{\link{bl_find_sparse_cf}}).
#'
#' @param bl_local_result A \code{"bl_local_result"} object produced by
#'   \code{\link{bl_find_local_cf}}.
#' @param exact_max_vars  Integer; use exact computation when the number of
#'   variables is at most this value. Default \code{14L}.
#' @param approx_perm     Integer; number of permutations for the approximate
#'   method. Default \code{2048L}.
#' @param seed            Random seed for the approximate method. \code{NULL}
#'   (default) uses the current RNG state, giving independent estimates across
#'   repeated calls. Supply an integer (e.g. \code{1L}) for reproducible results;
#'   this calls \code{set.seed()}, which resets the session's global RNG. Used
#'   only when the number of variables exceeds \code{exact_max_vars}; the exact
#'   method is deterministic and ignores \code{seed}.
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
#' @importFrom utils combn
#' @export
bl_shapley <- function(bl_local_result, exact_max_vars = 14L,
                       approx_perm = 2048L, seed = NULL) {
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
                               var_names)
  } else {
    shap <- .shapley_perm_one(start, end, bl_result$model, bl_result$model_type,
                              var_names,
                              M = as.integer(approx_perm), seed = seed)
  }
  names(shap) <- var_names

  pred_prob     <- bl_target$pred_prob
  pred_class    <- bl_target$pred_class
  pred_use      <- pred_class
  pred_boundary <- bl_local_result$B_pred
  data_to_boundary <- end - start

  # Ordering: class 0 -> decreasing (most negative first); class 1 -> ascending
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
    abs(df$shapley_cause) < 1e-8              ~ "No effect",
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
    round(df$pred_data, 3), " -> ", round(df$data_to_boundary, 2)
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
#' \emph{contradict} the prediction direction, or have \emph{no effect}
#' (dark blue, light grey and mid grey).
#'
#' @param x   A \code{"bl_shapley"} object from \code{\link{bl_shapley}}.
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot} object, returned visibly: it draws when called at
#'   the console. Inside a function, a loop or a \code{source()}d script, wrap
#'   the call in \code{print()} to draw it.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_label scale_fill_manual
#' @importFrom ggplot2 scale_color_manual theme_bw theme labs
#' @export
plot.bl_shapley <- function(x, ...) {
  df        <- x$shapley_df
  row_label <- if (is.na(x$row_id)) "external" else as.character(x$row_id)

  # Raw-unit relabel (display only): when a scaling was recorded, show each
  # variable's observed and counterfactual value in original units instead of
  # the default "observed -> change" in standardised units. Levels convert as
  # value*scale + center; the counterfactual level is pred_data + data_to_boundary.
  scaling <- x$bl_local_result$bl_result$scaling
  y_title <- "Variable (observed value -> change required)"
  if (!is.null(scaling)) {
    vn      <- as.character(df$varnames)
    obs_raw <- .scale_to_raw(df$pred_data, scaling, vn, "level")
    cf_raw  <- .scale_to_raw(df$pred_data + df$data_to_boundary,
                             scaling, vn, "level")
    lbl     <- paste0(vn, ": ", round(obs_raw, 3), " -> ", round(cf_raw, 2))
    df$varnames_p <- factor(lbl, levels = lbl)   # preserve existing sort order
    y_title <- "Variable (observed -> counterfactual, raw units)"
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(fill = Contribute, colour = Contribute,
                 x = shapley_cause, y = varnames_p)
  ) +
    ggplot2::scale_fill_manual(
      name   = "Contribution",
      values = c("Supports" = "darkblue", "Contradicts" = "grey70",
                 "No effect" = "grey45", "Unknown"  = "lightgrey")
    ) +
    ggplot2::scale_color_manual(
      name   = "Contribution",
      values = c("Supports" = "darkblue", "Contradicts" = "grey70",
                 "No effect" = "grey45", "Unknown"  = "lightgrey")
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
      y        = y_title,
      title    = "Shapley Contribution Plot",
      subtitle = sprintf(
        "ID: %s | Class: %d | Pred: %.3f | Boundary pred: %.3f",
        row_label, x$pred_class, x$pred_prob, x$pred_boundary
      )
    )
}


#' Print method for bl_shapley
#'
#' Prints the predictions at the observed point and the counterfactual, the sum
#' of the Shapley values, then a table in the same layout as \code{print()} of
#' a \code{bl_local_result}: one row per feature with the observed and
#' counterfactual values, the change, the \code{\link{set_filters}} constraint,
#' the Shapley value and its \code{Contribute} class. Feature values are in
#' original (raw) units when a scaling was set via \code{\link{bl_set_scaling}},
#' otherwise in model units; Shapley values are in prediction units. Numbers
#' print in fixed notation with 2 decimals, with more (up to 4) where a change
#' or a Shapley value would otherwise print as 0.00.
#'
#' @param x A \code{"bl_shapley"} object.
#' @param ... Currently ignored.
#' @export
print.bl_shapley <- function(x, ...) {
  cat("-- bl_shapley --\n")
  rid <- if (is.na(x$row_id)) "external" else as.character(x$row_id)
  cat(sprintf("  Row ID         : %s\n",  rid))
  cat(sprintf("  Observed pred  : %.4f  (class %d)\n", x$pred_prob, x$pred_class))
  cat(sprintf("  Full CF pred   : %.4f\n", x$pred_boundary))
  cat(sprintf("  Shapley sum    : %.4f\n", sum(x$shapley_df$shapley_cause)))

  df        <- x$shapley_df
  bl_local  <- x$bl_local_result
  var_names <- bl_local$bl_result$var_names
  x_obs     <- as.numeric(unlist(bl_local$bl_target$x_obs[1L, var_names]))
  B_x       <- as.numeric(unlist(bl_local$B_x[1L, var_names]))

  cmp <- .cf_compare_table(x_obs, list(counterfactual = B_x), var_names,
                           bl_local$bl_result$scaling,
                           bl_local$set_filters$constraints)
  idx <- match(var_names, as.character(df$varnames))
  s   <- df$shapley_cause[idx]
  tab <- .format_cf_table(cmp)
  tab$shapley    <- .fmt_fixed(s, .fmt_digits(s))
  tab$contribute <- df$Contribute[idx]
  cat(sprintf("\n  Observed vs counterfactual, with Shapley values (%s units):\n",
              cmp$units))
  print(tab)
  invisible(x)
}


#' Create a sparse counterfactual by zeroing contradicting variables
#'
#' Builds a sparse version of the counterfactual that changes as few variables
#' as needed. Starting from the observed point, it changes the variables whose
#' Shapley contribution \emph{supports} the prediction direction one at a time,
#' largest absolute Shapley value first, and re-scores the model after each.
#' It stops at the first point whose predicted class differs from the
#' target's; the remaining variables stay at their observed values. Variables
#' classified as \code{"Contradicts"}, \code{"No effect"} or \code{"Unknown"},
#' and any variable constrained as \code{"fixed"} in \code{\link{set_filters}},
#' are never changed. If no step crosses the boundary, every supporting
#' variable has been changed and \code{solution_valid} is \code{FALSE}.
#'
#' The selection follows the Shapley order, so it is not guaranteed to be the
#' smallest possible subset: a smaller set outside that order may also cross.
#'
#' @param bl_shapley_result A \code{"bl_shapley"} object from
#'   \code{\link{bl_shapley}}.
#' @param round_to Positive numeric or \code{NULL}. If not \code{NULL}, rounds
#'   each changed variable's value to the nearest multiple of this value (e.g.
#'   \code{round_to = 0.5} rounds to the nearest 0.5; \code{round_to = 1}
#'   rounds to the nearest integer). Rounding acts on the stored model-unit
#'   values, which are z-scores for a standardised model, and constraints are
#'   not re-checked afterwards; see
#'   \file{documentation/sparse-counterfactual.md}. Default \code{NULL} (no
#'   rounding).
#'
#' @return A list of class \code{"bl_sparse_result"} with the following fields:
#'   \describe{
#'     \item{x_sparse}{Single-row data frame; the sparse counterfactual in
#'       X-space.}
#'     \item{pred_sparse}{Scalar predicted probability of the sparse CF.}
#'     \item{solution_valid}{Logical; \code{TRUE} if the sparse CF flips the
#'       predicted class.}
#'     \item{shapley_df}{The annotated Shapley data frame with columns
#'       \code{used_in_sparse} (\code{TRUE} when the sparse CF changes the
#'       variable to its counterfactual value) and \code{x_sparse_val}.}
#'     \item{sparse_path}{Data frame, one row per step: \code{step},
#'       \code{variable} (the variable changed), \code{shapley} (its Shapley
#'       value), \code{pred} (the model probability after the step) and
#'       \code{crosses} (whether the predicted class has flipped). Zero rows
#'       when no variable is a candidate.}
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

  # "fixed"-constrained variables never change in the sparse CF, regardless
  # of Shapley classification (set_filters stored in bl_local_result).
  filters    <- bl_local$set_filters
  fixed_vars <- character(0)
  if (!is.null(filters) && length(filters$constraints) > 0L) {
    for (nm in names(filters$constraints)) {
      v <- filters$constraints[[nm]]
      if (is.character(v) && v == "fixed") fixed_vars <- c(fixed_vars, nm)
    }
  }

  # Candidates: "Supports" variables that are not fixed, largest |Shapley|
  # first ("Contradicts", "No effect" and "Unknown" stay at observed).
  vn_df  <- as.character(df$varnames)
  is_sup <- df$Contribute == "Supports" & !(vn_df %in% fixed_vars)
  sup    <- df[is_sup, , drop = FALSE]
  sup    <- sup[order(abs(sup$shapley_cause), decreasing = TRUE), , drop = FALSE]

  score <- function(x) {
    x_df <- as.data.frame(t(x))
    colnames(x_df) <- var_names
    .pred_function(
      model_use  = bl_result$model,
      model_type = bl_result$model_type,
      new_data   = x_df
    )[[1L]]
  }
  pred_class_target <- bl_target$pred_class

  # Greedy path: change one candidate at a time and re-score; stop at the
  # first point whose predicted class differs from the target's.
  x_sparse    <- x_obs
  pred_sparse <- score(x_sparse)
  path_rows   <- list()
  for (k in seq_len(nrow(sup))) {
    nm      <- as.character(sup$varnames[k])
    new_val <- B_x[[nm]]
    if (!is.null(round_to)) new_val <- round(new_val / round_to) * round_to
    x_sparse[[nm]] <- new_val
    pred_sparse    <- score(x_sparse)
    crosses        <- as.integer(pred_sparse >= bl_result$cutoff) != pred_class_target
    path_rows[[k]] <- data.frame(step = k, variable = nm,
                                 shapley = sup$shapley_cause[k],
                                 pred = pred_sparse, crosses = crosses,
                                 stringsAsFactors = FALSE)
    if (crosses) break
  }
  sparse_path <- if (length(path_rows) > 0L) {
    do.call(rbind, path_rows)
  } else {
    data.frame(step = integer(0), variable = character(0), shapley = numeric(0),
               pred = numeric(0), crosses = logical(0), stringsAsFactors = FALSE)
  }

  x_sparse_df <- as.data.frame(t(x_sparse))
  colnames(x_sparse_df) <- var_names

  pred_class_sparse <- as.integer(pred_sparse >= bl_result$cutoff)
  solution_valid    <- (pred_class_sparse != pred_class_target)

  # Annotate data frame
  df_out                <- df
  df_out$used_in_sparse <- vn_df %in% sparse_path$variable
  df_out$x_sparse_val   <- x_sparse[vn_df]

  structure(
    list(
      x_sparse       = x_sparse_df,
      pred_sparse    = pred_sparse,
      solution_valid = solution_valid,
      shapley_df     = df_out,
      sparse_path    = sparse_path,
      bl_shapley     = bl_shapley_result
    ),
    class = "bl_sparse_result"
  )
}


#' Plot a sparse counterfactual result
#'
#' Renders the rotated local biplot via \code{\link{plot.bl_local_result}}
#' (which uses the biplotEZ pipeline) and then overlays the sparse
#' counterfactual as a second cross marker. A green (\code{"green3"}) cross
#' indicates the sparse CF successfully flips the predicted class; a yellow
#' (\code{"yellow2"}) cross indicates it does not.
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
#'   \code{label_offset_var}, \code{which}, \code{new_title}).
#'   \code{rotate_deg} is not accepted: like the local biplot, the sparse
#'   counterfactual plot keeps the orientation set by the optimal rotation.
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
#' Prints the predictions at the observed point, the full counterfactual and
#' the sparse counterfactual, then a table in the same layout as
#' \code{print()} of a \code{bl_local_result}: one row per feature with the
#' observed, counterfactual and sparse values, the change the sparse
#' counterfactual asks for (sparse minus observed), the
#' \code{\link{set_filters}} constraint, and whether the sparse counterfactual
#' uses the feature's counterfactual value. Values are in original (raw) units
#' when a scaling was set via \code{\link{bl_set_scaling}}, otherwise in model
#' units. Last, it prints the path \code{\link{bl_find_sparse_cf}} took: the
#' feature changed at each step, its Shapley value, the prediction after the
#' step, and whether the predicted class has flipped. Numbers print in fixed
#' notation with 2 decimals, with more (up to 4) where a change or Shapley value
#' would otherwise print as 0.00, or where a prediction is within 0.005 of the
#' cutoff.
#'
#' @param x A \code{"bl_sparse_result"} object.
#' @param ... Currently ignored.
#' @export
print.bl_sparse_result <- function(x, ...) {
  cat("-- bl_sparse_result --\n")
  cat(sprintf("  Solution valid : %s\n", x$solution_valid))
  cat(sprintf("  Observed pred  : %.4f  (class %d)\n",
              x$bl_shapley$pred_prob, x$bl_shapley$pred_class))
  cat(sprintf("  Full CF pred   : %.4f\n", x$bl_shapley$pred_boundary))
  cat(sprintf("  Sparse CF pred : %.4f\n", x$pred_sparse))

  df        <- x$shapley_df
  bl_local  <- x$bl_shapley$bl_local_result
  var_names <- bl_local$bl_result$var_names
  x_obs     <- as.numeric(unlist(bl_local$bl_target$x_obs[1L, var_names]))
  B_x       <- as.numeric(unlist(bl_local$B_x[1L, var_names]))
  x_sparse  <- as.numeric(unlist(x$x_sparse[1L, var_names]))

  cmp <- .cf_compare_table(x_obs, list(counterfactual = B_x, sparse = x_sparse),
                           var_names, bl_local$bl_result$scaling,
                           bl_local$set_filters$constraints)
  tab      <- .format_cf_table(cmp)
  tab$used <- df$used_in_sparse[match(var_names, as.character(df$varnames))]
  cat(sprintf("\n  Observed vs sparse counterfactual (%s units):\n", cmp$units))
  print(tab)
  cat("  (change = sparse - observed)\n")

  path <- x$sparse_path
  if (!is.null(path)) {
    if (nrow(path) == 0L) {
      cat("\n  Path: no \"Supports\" feature to change.\n")
    } else {
      cat("\n  Path (Supports features, largest Shapley value first):\n")
      # A prediction within 0.005 of the cutoff keeps enough decimals to
      # show which side it is on.
      cutoff       <- bl_local$bl_result$cutoff
      path$shapley <- .fmt_fixed(path$shapley, .fmt_digits(path$shapley))
      path$pred    <- .fmt_fixed(path$pred, .fmt_digits(path$pred - cutoff))
      print(path, row.names = FALSE)
    }
  }
  invisible(x)
}
