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
#'   relative to within-class spread. CVA technically requires at least 3
#'   distinct class levels, but biplotEZ supports 2-class data (results may
#'   be degenerate). When `bl_model` is supplied, the four confusion
#'   categories (TP, TN, FP, FN) are used as the class labels, satisfying
#'   the 3-class requirement. As a simpler alternative, pass
#'   `cva_classes = as.factor(train_data$class)` directly (2 levels only).
#'
#' @param train_data  Data frame with feature columns and a `"class"` column
#'   (numeric 0/1). The class column is used to derive confusion labels for
#'   CVA when `bl_model` is supplied.
#' @param var_names   Character vector of feature column names.
#' @param method      `"PCA"` or `"CVA"`. Default `"PCA"`.
#' @param standardise Logical; whether to standardise features before PCA.
#'   Ignored (forced to `FALSE`) when `method = "CVA"`. Default `TRUE`.
#' @param proj_dims   Integer vector of length 2; which eigenvector indices
#'   to use as the biplot plane. Default `c(1L, 2L)`.
#' @param bl_model    A `"bl_model"` object from `bl_fit_model()` or
#'   `bl_wrap_model()`. When `method = "CVA"` and `cva_classes` is `NULL`,
#'   `bl_model` is used to compute four-class confusion labels (TP, TN, FP,
#'   FN) which are passed as `cva_classes`. Ignored for PCA.
#' @param cva_classes Factor of class labels, one per row of `train_data`.
#'   When supplied, takes precedence over `bl_model` for CVA. Use
#'   `as.factor(train_data$class)` as a simple alternative when no model
#'   has been fitted.
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
#' bl_dat  <- bl_prepare_data(datasets::iris,
#'                             class_col    = "Species",
#'                             target_class = "versicolor")
#' # PCA projection
#' bl_proj <- bl_build_projection(bl_dat$train_data,
#'                                 bl_dat$var_names,
#'                                 method = "PCA")
#'
#' # CVA with bl_model — confusion labels (TP/TN/FP/FN) used automatically
#' bl_mod <- bl_fit_model(bl_dat$train_data, bl_dat$var_names)
#' bl_proj_cva <- bl_build_projection(
#'   bl_dat$train_data,
#'   bl_dat$var_names,
#'   method    = "CVA",
#'   proj_dims = c(1L, 2L),
#'   bl_model  = bl_mod
#' )
#'
#' # CVA alternative: pass binary class factor directly (2 levels only)
#' bl_proj_cva2 <- bl_build_projection(
#'   bl_dat$train_data,
#'   bl_dat$var_names,
#'   method      = "CVA",
#'   proj_dims   = c(1L, 2L),
#'   cva_classes = as.factor(bl_dat$train_data$class)
#' )
#'
#' @importFrom biplotEZ biplot PCA CVA
#' @importFrom dplyr case_when
#' @importFrom stats sd
#' @export
bl_build_projection <- function(train_data,
                                var_names,
                                method      = "PCA",
                                standardise = TRUE,
                                proj_dims   = c(1L, 2L),
                                bl_model    = NULL,
                                cva_classes = NULL,
                                title       = "") {

  # ---- Validation -------------------------------------------------------
  stop_if_not_data_frame(train_data, "train_data")
  stop_if_not_character(var_names, "var_names")
  method <- match.arg(method, c("PCA", "CVA"))

  if (!is.numeric(proj_dims) || length(proj_dims) != 2L)
    stop("'proj_dims' must be an integer vector of length 2.", call. = FALSE)
  proj_dims <- as.integer(proj_dims)

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
    # ---- CVA: derive cva_classes if not supplied -----------------------
    if (is.null(cva_classes)) {
      if (!is.null(bl_model) && inherits(bl_model, "bl_model")) {
        # Compute four-class confusion labels from the fitted model
        if (!"class" %in% names(train_data))
          stop(paste0(
            "method = 'CVA' with bl_model requires a 'class' column in ",
            "train_data (numeric 0/1) to compute confusion labels."
          ), call. = FALSE)

        
    ##A This is a custom prediction function that can be replaced by a standard R model    
        pred_prob  <- .pred_function(
          model_use  = bl_model$model,
          model_type = bl_model$model_type,
          new_data   = train_data[, var_names, drop = FALSE]
        )
        pred_class <- as.numeric(pred_prob >= bl_model$cutoff)
        actual     <- train_data[["class"]]

    ##A key output used to colour points and in CVA create the classes: just used for the data points
        confusion_label <- dplyr::case_when(
          actual == pred_class & actual == 1L ~ "TP",
          actual == pred_class & actual == 0L ~ "TN",
          actual != pred_class & actual == 0L ~ "FP",
          TRUE                               ~ "FN"
        )
        cva_classes <- as.factor(confusion_label)

      } else {
        stop(paste0(
          "method = 'CVA' requires either:\n",
          "  (a) bl_model: a 'bl_model' object — confusion labels (TP/TN/FP/FN) ",
          "will be used automatically, or\n",
          "  (b) cva_classes: a factor of class labels (one per row of train_data),\n",
          "      e.g., as.factor(train_data$class) when no model is available."
        ), call. = FALSE)
      }
    }

    if (length(levels(cva_classes)) < 3L)
      warning(paste0(
        "CVA works best with at least 3 class levels. ",
        "Only ", length(levels(cva_classes)), " level(s) found. ",
        "Consider supplying bl_model to use TP/TN/FP/FN labels (4 levels), ",
        "or use method = 'PCA'."
      ), call. = FALSE)

    bp <- biplotEZ::biplot(X, scaled = FALSE, Title = title) |>
      biplotEZ::CVA(e.vects          = proj_dims,
                    classes          = cva_classes,
                    show.class.means = FALSE)
  }

  # ---- Extract V and tV ------------------------------------------------
  V  <- bp$Lmat
  tV <- solve(V)

  # ---- Default point colours -------------------------------------------
  # Stored as metadata so plot.bl_projection() can colour automatically.
  # Computation only; no plotting happens here.
  # Per-observation colour vector. Used by plot.bl_projection() via
  # graphics::points(), which accepts per-obs colours directly (unlike
  # biplotEZ::samples(col=...) which maps colours to class levels).
  point_col <- if (method == "CVA" && !is.null(cva_classes)) {
    if (setequal(levels(cva_classes), c("FN", "FP", "TN", "TP"))) {
      col_map <- c(TP = "green3", TN = "steelblue", FP = "orange", FN = "red")
      unname(col_map[as.character(cva_classes)])
    } else {
      ifelse(as.character(cva_classes) == "1", "red", "blue")
    }
  } else if ("class" %in% names(train_data)) {
    ifelse(train_data[["class"]] == 1L, "red", "blue")
  } else {
    NULL
  }

  # ---- Training variable ranges (min/max per feature) ------------------
  train_ranges <- get_variable_ranges(train_data[, var_names, drop = FALSE])

  # ---- Return ----------------------------------------------------------
  structure(
    list(
      V            = V,
      tV           = tV,
      X_center     = X_center,
      X_sd         = X_sd,
      method       = method,
      standardise  = standardise,
      proj_dims    = proj_dims,
      biplot_obj   = bp,
      cva_classes  = cva_classes,
      point_col    = point_col,
      train_ranges = train_ranges
    ),
    class = "bl_projection"
  )
}


# --------------------------------------------------------------------------
# S3 plot method
# --------------------------------------------------------------------------

#' Plot a bl_projection biplot
#'
#' Renders the biplotEZ object stored in a `"bl_projection"` result.
#' Default point colours are derived automatically:
#' - CVA with TP/TN/FP/FN labels: green/steelblue/orange/red confusion colours.
#' - CVA with binary 0/1 labels: red (class 1) / blue (class 0).
#' - PCA (or any other case): red (class 1) / blue (class 0) when a
#'   `"class"` column was present in `train_data`; grey otherwise.
#'
#' @param x   A `"bl_projection"` object from [bl_build_projection()].
#' @param col Character vector of point colours (one per training row).
#'   When `NULL` (default) the colours stored in `x$point_col` are used.
#' @param pch Integer point character. Default `16`.
#' @param cex Numeric point size. Default `0.8`.
#' @param axes_col Colour for axis lines and labels. Default `"grey22"`.
#' @param tick_cex Numeric tick-label size. Default `0.6`.
#' @param ... Ignored; present for S3 compatibility.
#'
#' @return Invisibly returns `x`.
#' @export
plot.bl_projection <- function(x,
                                col      = NULL,
                                pch      = 16,
                                cex      = 0.8,
                                axes_col = "grey22",
                                tick_cex = 0.6,
                                ...) {
  col <- if (!is.null(col)) col else x$point_col
  if (is.null(col)) col <- "grey40"

  # Use the same pattern as plot_biplotEZ(): hide biplotEZ's own sample
  # rendering (opacity = 0) then draw points manually via graphics::points().
  # This avoids biplotEZ::samples(col=...) mapping colours to class levels
  # rather than to individual observations.
  x$biplot_obj |>
    biplotEZ::samples(opacity = 0, which = NULL) |>
    biplotEZ::axes(col            = axes_col,
                   label.dir      = "Hor",
                   tick.label.cex = tick_cex,
                   ticks          = 1L) |>
    plot()

  graphics::points(x$biplot_obj$Z[, x$proj_dims],
                   col = col,
                   pch = pch,
                   cex = cex)

  invisible(x)
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
