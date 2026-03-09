
############################################################
# Script: Able to select points in the biplot to identify data poin and associated CF
# Author: Adriaan Rowan
# Date: 2025-10-01
# Description: Pick points on the biplot and return row + data

#########################################################
#' Interactively pick points from the current biplot
#'
#' Works with base graphics. You must have already drawn the biplot points
#' (e.g., via point_inter_proj(..., plot = TRUE)) so the points are visible.
#'
#' @param X           The original data (same rows as projected).
#' @param proj        Projection indices (default pc.values).
#' @param standardise Logical; whether to standardize by Xsd.
#' @param CVA_method  Logical; if TRUE, force standardise = FALSE.
#' @param Xcenter     Numeric means used in projection.
#' @param Xsd         Numeric sds used in projection.
#' @param V           Loadings/eigenvector matrix.
#' @param n           How many points to pick (default 1). ESC to stop early.
#' @param method      "identify" (fast; nearest plotted symbol) or "locator"
#'                    (click anywhere; nearest point in 2D).
#' @param show_labels Whether to label picked points on the plot.
#' @param label_cex   Label size.
#'
#' @return A list with:
#'   - idx: integer vector of selected row indices
#'   - picked: data.frame with row index, projected coords, and original X row
pick_biplot_points <- function(
    X,
    proj = pc.values,
    standardise = FALSE, CVA_method = FALSE,
    Xcenter = X_center, Xsd = X_sd,
    V = V,
    n = 1,
    method = c("identify", "locator"),
    show_labels = TRUE,
    label_cex = 0.9
) {
  method <- match.arg(method)
  
  # Compute full projected coordinates (no subsetting)
  sv <- Xsd
  if (isTRUE(CVA_method)) standardise <- FALSE
  if (!isTRUE(standardise)) sv <- FALSE
  X.st <- scale(X, center = Xcenter, scale = sv)
  
  Z_all <- X.st %*% V[, proj, drop = FALSE]
  colnames(Z_all) <- c("x", "y")
  
  # Helper: assemble return data for chosen indices
  .assemble <- function(idxs) {
    idxs <- unique(stats::na.omit(as.integer(idxs)))
    if (length(idxs) == 0L) {
      return(list(idx = integer(0), picked = data.frame()))
    }
    df <- cbind(
      data.frame(row = idxs,
                 z_x = Z_all[idxs, 1],
                 z_y = Z_all[idxs, 2]),
      as.data.frame(X[idxs, , drop = FALSE])
    )
    rownames(df) <- NULL
    if (isTRUE(show_labels)) {
      text(Z_all[idxs, 1], Z_all[idxs, 2], labels = idxs, pos = 3, cex = label_cex)
    }
    list(idx = idxs, picked = df)
  }
  
  if (method == "identify") {
    # Requires points already drawn on the active device
    picked_idx <- graphics::identify(x = Z_all[, 1], y = Z_all[, 2],
                                     labels = seq_len(nrow(Z_all)),
                                     n = n, plot = FALSE)
    return(.assemble(picked_idx))
  } else {
    # "locator": click anywhere; choose nearest point each time
    idxs <- integer(0)
    for (k in seq_len(n)) {
      loc <- tryCatch(graphics::locator(1), error = function(e) NULL)
      if (is.null(loc)) break
      d2 <- (Z_all[, 1] - loc$x)^2 + (Z_all[, 2] - loc$y)^2
      idxs <- c(idxs, which.min(d2))
    }
    return(.assemble(idxs))
  }
}
