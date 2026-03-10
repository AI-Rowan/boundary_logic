############################################################
# bl_build_projection(): PCA or CVA biplot projection matrix
# Refactored from: scripts/1 Functions for biplot.R (biplot_input_calc)
############################################################

#' Build a PCA or CVA projection matrix for boundary logic analysis
#'
#' Computes the loading matrix V, its inverse tV, and centering/scaling
#' parameters from the training data. Uses `biplotEZ` for the PCA/CVA
#' computation and biplot object construction.
#'
#' @section CVA standardisation:
#' CVA always forces `standardise = FALSE` internally, regardless of the
#' value supplied by the user. This is because CVA uses within-class
#' covariance structure that already accounts for scale. The effective
#' `standardise` value is recorded in the returned object.
#'
#' @section Choosing between PCA and CVA:
#' - **PCA** (`method = "PCA"`): suitable for any data; separates
#'   observations by overall variance. Use when class structure is not
#'   well-defined in the feature space.
#' - **CVA** (`method = "CVA"`): maximises between-class separation
#'   relative to within-class spread. Requires at least 3 distinct class
#'   labels in `cva_classes` for a 2D projection.
#'
#' @param train_data  Data frame; feature columns only (no class column).
#' @param var_names   Character vector of feature column names.
#' @param method      `"PCA"` or `"CVA"`. Default `"PCA"`.
#' @param standardise Logical; whether to standardise features before PCA.
#'   Ignored (forced to `FALSE`) when `method = "CVA"`. Default `TRUE`.
#' @param proj_dims   Integer vector of length 2; which eigenvector indices
#'   to use as the biplot plane. Default `c(1L, 2L)`.
#' @param cva_classes Factor of class labels, one per row of `train_data`.
#'   Required when `method = "CVA"`; ignored for PCA.
#' @param title       Character string for the biplot title. Default `""`.
#'
#' @return A list of class `"bl_projection"` with components:
#' \describe{
#'   \item{`V`}{Numeric matrix (p x p); loading / eigenvector matrix.}
#'   \item{`tV`}{Numeric matrix (p x p); `solve(V)`.}
#'   \item{`X_center`}{Numeric vector (length p); column means of the
#'     feature matrix.}
#'   \item{`X_sd`}{Numeric vector (length p); column standard deviations.}
#'   \item{`method`}{`"PCA"` or `"CVA"`.}
#'   \item{`standardise`}{Logical; the effective standardise flag after any
#'     CVA override.}
#'   \item{`proj_dims`}{Integer vector; the eigenvector indices used.}
#'   \item{`biplot_obj`}{The `biplotEZ` S3 object; used for plotting in
#'     Phase 2.}
#' }
#'
#' @examples
#' bl_dat  <- bl_prepare_data(iris,
#'                             class_col    = "Species",
#'                             target_class = "versicolor")
#' # PCA projection
#' bl_proj <- bl_build_projection(bl_dat$train_data,
#'                                 bl_dat$var_names,
#'                                 method = "PCA")
#'
#' # CVA projection (requires class labels with >= 3 levels for 2D)
#' pred_col <- ifelse(bl_dat$train_data$class == 1, "red", "blue")
#' bl_proj_cva <- bl_build_projection(
#'   bl_dat$train_data,
#'   bl_dat$var_names,
#'   method      = "CVA",
#'   cva_classes = as.factor(pred_col)
#' )
#'
#' @export
bl_build_projection <- function(train_data,
                                var_names,
                                method      = "PCA",
                                standardise = TRUE,
                                proj_dims   = c(1L, 2L),
                                cva_classes = NULL,
                                title       = "") {

  # ---- Validation -------------------------------------------------------
  stop_if_not_data_frame(train_data, "train_data")
  stop_if_not_character(var_names, "var_names")
  method <- match.arg(method, c("PCA", "CVA"))

  if (!is.numeric(proj_dims) || length(proj_dims) != 2L)
    stop("'proj_dims' must be an integer vector of length 2.", call. = FALSE)
  proj_dims <- as.integer(proj_dims)

  if (method == "CVA" && is.null(cva_classes))
    stop(paste0(
      "method = 'CVA' requires 'cva_classes' to be supplied. ",
      "Pass a factor of class labels (one per row of train_data), ",
      "e.g., as.factor(train_data$class) or the prediction colour vector."
    ), call. = FALSE)

  # ---- CVA forces standardise = FALSE ----------------------------------
  if (method == "CVA" && isTRUE(standardise)) {
    standardise <- FALSE
  }

  # ---- Feature matrix --------------------------------------------------
  X        <- as.matrix(train_data[, var_names, drop = FALSE])
  X_center <- colMeans(X)
  X_sd     <- apply(X, 2L, stats::sd)

  # ---- Build biplotEZ object -------------------------------------------
  if (method == "PCA") {
    bp <- biplotEZ::biplot(X, scaled = standardise, Title = title) |>
      biplotEZ::PCA(e.vects = proj_dims)
  } else {
    # CVA: biplotEZ expects the original data frame, not scaled
    if (length(levels(cva_classes)) < 3L)
      warning(paste0(
        "CVA requires at least 3 class levels for a 2D projection. ",
        "Only ", length(levels(cva_classes)), " level(s) found. ",
        "Results may be degenerate. Consider using method = 'PCA'."
      ), call. = FALSE)

    bp <- biplotEZ::biplot(X, scaled = FALSE, Title = title) |>
      biplotEZ::CVA(e.vects      = proj_dims,
                    classes      = cva_classes,
                    show.class.means = FALSE)
  }

  # ---- Extract V and tV ------------------------------------------------
  V  <- bp$Lmat
  tV <- solve(V)

  # ---- Return ----------------------------------------------------------
  structure(
    list(
      V           = V,
      tV          = tV,
      X_center    = X_center,
      X_sd        = X_sd,
      method      = method,
      standardise = standardise,
      proj_dims   = proj_dims,
      biplot_obj  = bp
    ),
    class = "bl_projection"
  )
}


# --------------------------------------------------------------------------
# S3 print method
# --------------------------------------------------------------------------

#' @export
print.bl_projection <- function(x, ...) {
  cat("<bl_projection>\n")
  cat(sprintf("  Method      : %s\n", x$method))
  std_note <- if (x$method == "CVA") {
    "FALSE (forced by CVA)"
  } else {
    as.character(x$standardise)
  }
  cat(sprintf("  Standardise : %s\n", std_note))
  cat(sprintf("  Dimensions  : %s\n", paste(x$proj_dims, collapse = ", ")))
  cat(sprintf("  V matrix    : %d x %d\n", nrow(x$V), ncol(x$V)))
  invisible(x)
}
