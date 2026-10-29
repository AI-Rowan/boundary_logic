# Originated in the `Explain XAI code` project's sandbox as
# `R/In-development/17_2_loan-gam-pca-4var-slices.R` (snapshot 2026-09-10). This
# repo is now this script's home; edit it here.

# ============================================================================
# Why two rotations give different decision boundaries — loan data, 4 variables
# ----------------------------------------------------------------------------
# WHY THIS SCRIPT EXISTS
#
# Script 17_1 makes the rotation confirmable by dropping to THREE features, where
# PCA gives a 3 x 3 orthogonal V and the 3D view is exactly lossless. Its central
# result (its section 18) is that a 2D biplot is a SLICE of the decision surface,
# not a view of it: two planes through one target cut one surface along two
# genuinely different curves, so a different plane gives a different
# counterfactual and different advice.
#
# This script carries that result to FOUR features, where it matters more and
# where three things are true that were not true at p = 3.
#
# WHAT CHANGED FROM THE PREVIOUS VERSION OF THIS SCRIPT (2026-09-04)
#
# The green GAM "decision surface" is GONE, along with the p ~ 0.5 cloud that fed
# it. It was fitted as gam(s3 ~ s(s1, s2, k = 60)) through a band of grid nodes
# with |p - 0.5| < 0.02 and then drawn in a 3-slot view of a 4-dimensional space.
# The previous header listed five reasons it disagreed with the boundary Boundary
# Logic actually computes -- smoothing band, single-valuedness in s3, sparse
# support at the edges, three separate discretisations, and a slot-4 pinning that
# differed from the package's. All five are properties of the FITTED SURFACE, not
# of the method. Removing it does not lose information; it removes a fabrication.
#
# What is left is exact. The package's contours ARE the p = 0.5 curves of the real
# 4D GAM restricted to the drawn plane -- contoured on a dense grid, not smoothed
# through a cloud. So the boundary is drawn only where it is exact:
#     view A -> bl_res$biplot_grid$ct   (the global biplot's boundary, in PC1-PC2)
#     view B -> loc_12$ct_local         (pinned to pair (1,2))
#     view C -> loc_op$ct_local         (the pair the package selects)
# all plotted verbatim, in the plane's own 2D coordinates.
#
# ----------------------------------------------------------------------------
# THE ENABLING RESULT: A 3D FIGURE CAN BE EXACT HERE
#
# Both rotated planes contain the target, and they are distinct. Two 2-planes
# sharing exactly a line span 2 + 2 - 1 = 3 dimensions. So there is a
# 3-DIMENSIONAL SUBSPACE CONTAINING BOTH PLANES EXACTLY, and a 3D figure drawn in
# it holds both planes, both boundary curves, the target and both counterfactuals
# with ZERO approximation. Only the 221 data points are projected, and that loss
# is measured in section 17.
#
# Section 12 builds that subspace as W = [Vr_B, n3], with n3 the single direction
# plane C adds to plane B. Choosing the frame this way makes plane B the z = 0
# plane, so plane C is seen tilting about the shared ray. The proof that the two
# planes span exactly three dimensions is one line: the residual of plane C after
# removing plane B has singular values (4.54e-01, 1.42e-15) -- rank ONE.
#
# ----------------------------------------------------------------------------
# THREE THINGS THAT ARE TRUE AT p = 4 AND WERE NOT AT p = 3
#
# 1. THE SHARED RAY IS A REAL RESULT, NOT A DIMENSIONAL ACCIDENT.
#    In 3-space any two 2-planes must intersect in a line, so 17_1's shared ray
#    was forced by dimension and carried no information -- its section 12 says so
#    explicitly, and its first principal angle was 0 for every pair for that
#    reason. In 4-space two 2-planes generically meet only at the ORIGIN. That
#    these two meet in a line, and that the line is the target's own ray, is
#    therefore genuine content. Section 12 computes the intersection WITHOUT
#    assuming the answer (the null space of [Vr_B, -Vr_C]) and then checks it:
#    exactly one zero singular value, and |cos(intersection, x_t)| = 1 to 12 dp.
#    This is the p = 4 home of the reference note's section 4.1 remark that the
#    first principal angle marks the target's direction.
#
# 2. THE 3D PICTURE CAN MISLEAD, AND SECTION 16 MEASURES HOW.
#    At p = 3 every basis was a lossless rotation, so no 3D view could mislead.
#    Here the unrotated PC1-PC2-PC3 view drops PC4, and PC4 is a direction the
#    planes genuinely use: the curves carry a PC4 component of up to 2.93 against a
#    data span of 7.155.
#
#    BE PRECISE ABOUT WHAT THAT BREAKS -- a first pass at section 16 got this wrong
#    and the correction is worth carrying. Projection is LINEAR, so it preserves
#    incidence: planes stay planes, curves stay in their planes, and the point the
#    two curves must share stays shared (measured: its distance to each curve moves
#    from 0.0166 to 0.0160, i.e. not at all). The figure is not topologically
#    broken. What projection destroys is the METRIC. DISTANCES can only shrink, and
#    the worst-hit distance is the one the whole method is about: the target's
#    distance to its counterfactual -- how much would have to change for the
#    decision to flip -- reads 0.281 as drawn when it is 0.657, i.e. 57% short.
#    ANGLES are not distances and can move either way: the angle between the two
#    planes comes out slightly LARGER as drawn (28.02 vs 26.98 deg). So "the
#    projection shrinks everything" would be wrong; it distorts direction by
#    direction, and nothing measured off the PCA figure can be trusted in either
#    direction. Section 16 tabulates each quantity, true versus as-drawn, in both
#    bases.
#
#    The trade inverts the obvious ranking, which is why it is worth stating: the
#    PCA basis represents the DATA CLOUD better (median lost fraction 0.045 versus
#    0.182 in W -- PCA maximises retained variance, that is its job) and the
#    GEOMETRY worse. "Better" has to name what it is better AT.
#
# 3. THE UNDER-DETERMINATION IS WORSE, AND ITS COST IS MEASURED WITHIN ONE PAIR.
#    Section 8. The Procrustes matrix is rank 1 of 4, so 3 of the rotation's
#    degrees of freedom are an arbitrary completion (at p = 3 only 1 was). Section
#    8 measures the free dimensions rather than quoting a ladder -- see the CAREFUL
#    note there, because the count of free ROTATIONS and the count of free PLANES
#    are different numbers and the reference note's section 5.1 conflates them.
#
# ----------------------------------------------------------------------------
# HOW TO READ THE ROTATED VIEWS: THEY ARE **RC BIPLOTS**, NOT "ROTATED PCs"
#
# After the rotation the axes are RC1 / RC2 -- "rotated components", the thesis's
# own term for a rotated component (RC) biplot. They are NOT "rotated PC2 / PC3",
# and the named pair is NOT an axis identity. From the reference note
# (biplot_rotation_and_pair_selection.md, section 2.7):
#
#   "the SVD construction guarantees the target's rotated information concentrates
#    in columns 1-2 of Vrho REGARDLESS of which eigenvector pair was used to build
#    Vr -- which is why the function can simply always return columns 1-2 rather
#    than tracking where the pair's original indices ended up."
#
# The previous version of this script labelled its axes "Rotated PC%d", which is
# exactly what that paragraph forbids. Section 9 prints the evidence for whatever
# target is set: per pair, the rotation angle theta, the principal angles against
# the eigenplane the pair NAMES, how much of the rotated plane actually lies inside
# that eigenplane, and which original PC each rotated column is closest to. The
# named pair is provenance for WHICH PLANE was rotated, so it goes in figure
# titles, never on an axis.
#
# ----------------------------------------------------------------------------
# TERMINOLOGY DISCIPLINE
#
# "Optimal" and "best" refer ONLY to what bl_find_local_cf() measures: the
# smallest squared Mahalanobis distance from the target to its counterfactual,
# under the inverse pooled within-class covariance. How steeply a plane is tilted,
# how much variance a basis retains, and how readable a figure is are NOT criteria
# the package optimises, and they can point the opposite way (see finding 2 above,
# where the basis that keeps more variance is the one that breaks the geometry).
#
# EVERYTHING GEOMETRIC COMES FROM boundarylogic. The PCs are bl_res$V; the
# rotations are bl_find_local_cf()$Vr_rot; the boundaries are ct_local and
# biplot_grid$ct; cutoff, b_margin and metric_inv come from bl_res. svd(Xstd)$v
# appears nowhere in the drawn path -- section 7 explains why that matters here
# and is not a stylistic preference.
# ============================================================================

rm(list = ls())

library(boundarylogic)
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
if (!requireNamespace("mgcv",   quietly = TRUE)) install.packages("mgcv")
library(plotly)
library(mgcv)

# --- 0) Analysis options ----------------------------------------------------
# Target row. A usable target must satisfy three hard requirements:
#   1. the package's selected pair is not (1,2), or views B and C are the same
#      search and the comparison shows nothing;
#   2. the B and C PLANES are materially apart -- checked by the SECOND PRINCIPAL
#      ANGLE between them, NOT by the pair label. Measured on this data, rows 12
#      (3,4), 61 (1,2), 77 (1,4) and 113 (1,2) all carry a pair label while their
#      B and C planes coincide EXACTLY (0.00 deg), so the two views would be one
#      picture drawn twice. This is 17_1 section 11's collapse guard, and it fires
#      on the plane, not the label;
#   3. best_pair is reproducible run to run. The note's section 5.3 records row 102
#      flipping between (1,2) and (1,3) across five identical calls at p = 3. All
#      15 rows scanned here were stable at p = 4, but the check is kept and run.
#
# ROW 49 (set below) satisfies all three: selected pair (2,3), planes 26.98 deg
# apart, stable across five identical calls. It is also the row the reference
# note's section 4.1 perturbation table is measured on, so this script stays
# continuous with the note.
#
# Alternatives, if a more extreme contrast is wanted: rows 6 (1,3), 95 (1,4) and
# 175 (2,3) all give planes a full 90.00 deg apart -- orthogonal, the largest
# separation possible. Row 104 (2,4) and row 33 (2,4) give 45.00 deg.
i <- 49L

feat_names <- c("loan_percent_income", "credit_score", "loan_int_rate", "loan_amnt")

m_probe <- 200L   # grid resolution for section 8's measurement recipe (matches package)

# --- 1) Load + filter loan data ---------------------------------------------
loan_raw <- read.csv(
  system.file("extdata", "loan_data.csv", package = "boundarylogic")
)
loan_filtered <- loan_raw[loan_raw$previous_loan_defaults_on_file == "No", ]
loan_filtered <- loan_filtered[, c(feat_names, "loan_status")]
loan_filtered$credit_score <- as.numeric(loan_filtered$credit_score)
cat("Retained rows:", nrow(loan_filtered), "\n")

# --- 2) Prepare data --------------------------------------------------------
bl_dat <- bl_prepare_data(
  data           = loan_filtered,
  class_col      = "loan_status",
  feature_cols   = feat_names,
  train_fraction = 0.01,
  seed           = 121L,
  hull_fraction  = 1
)
data <- bl_dat$train_data

# --- 3) Standardize the 4 features ------------------------------------------
sc       <- scale(data[, feat_names], center = TRUE, scale = TRUE)
X_center <- attr(sc, "scaled:center")
X_scale  <- attr(sc, "scaled:scale")
Xstd     <- sc                              # n x 4, standardized
p        <- ncol(Xstd)
stopifnot(i >= 1L, i <= nrow(Xstd))
x_t <- as.numeric(Xstd[i, ])

# --- 4) Fit the nonlinear GAM in standardized 4D space ----------------------
dat_gam <- data.frame(Xstd, class = data$class)
names(dat_gam)[1:4] <- feat_names
gam_fit <- gam(
  class ~ s(loan_percent_income, k = 8) + s(credit_score, k = 8) +
          s(loan_int_rate, k = 8) + s(loan_amnt, k = 8),
  data = dat_gam, family = binomial()
)
summary(gam_fit)

pt_col <- ifelse(dat_gam$class == 1, "red", "blue")   # colour by ACTUAL class

# Model probability at a matrix of standardised feature-space points.
pfun <- function(X) {
  X <- as.data.frame(X); colnames(X) <- feat_names
  predict(gam_fit, newdata = X, type = "response")
}

# --- 5) The boundarylogic pipeline — the source of the PCs and the metric ----
# bl_fit_model() cannot fit a GAM -- GAM/GBM/LDA/XGB were removed from it (see the
# package NEWS.md); bl_wrap_model() is the supported route, and needs no
# predict_fn because "GAM" is natively recognised and predicted with
# type = "response".
#
# The data handed to bl_prepare_data() is standardised in place here, so the
# package's "raw" feature space IS this script's Xstd space. That makes every
# comparison below a direct one rather than a units conversion.
bl_dat$train_data[, feat_names] <- Xstd
bl_dat$test_data[,  feat_names] <-
  scale(bl_dat$test_data[, feat_names], center = X_center, scale = X_scale)
bl_dat <- bl_set_scaling(bl_dat, center = X_center, scale = X_scale,
                         method = "z-score")

bl_mod <- bl_wrap_model(gam_fit, "GAM", bl_dat$var_names,
                        train_data = bl_dat$train_data)
bl_res <- bl_build_result(bl_data = bl_dat, bl_model = bl_mod,
                          method = "PCA", title = "Loan Status (4 vars)")

V_full     <- bl_res$V           # THE PCs: 4 x 4, the package's loading matrix
cutoff     <- bl_res$cutoff      # 0.5
b_margin   <- bl_res$b_margin    # 0.001
metric_inv <- bl_res$metric_inv  # inverse pooled within-class covariance

stopifnot(max(abs(crossprod(V_full) - diag(p))) < 1e-10)
cat(sprintf("\nPackage V vs svd(Xstd)$v: max abs diff %.3e (section 7: this matters)\n",
            max(abs(V_full - svd(Xstd)$v))))
if (is.null(metric_inv))
  warning("bl_result carries no metric_inv -- pair selection falls back to Euclidean.",
          call. = FALSE)

# --- 6) The two planes, both from the package -------------------------------
# max_pairs = 1L pins the search to pair (1,2), because .generate_pairs() emits
# colexicographically -- for j in 2..p, for i in 1..j-1, so (1,2) is always first.
# It is the only supported way to force a plane, and it is what makes a
# pinned-versus-selected comparison possible at all.
tgt    <- bl_select_target(bl_res, target = bl_res$train_data[i, ])
loc_12 <- bl_find_local_cf(bl_res, tgt, max_pairs = 1L, verbose = FALSE)
loc_op <- bl_find_local_cf(bl_res, tgt, verbose = FALSE)

if (!isTRUE(loc_12$solution_found) || !isTRUE(loc_op$solution_found))
  stop("bl_find_local_cf() found no counterfactual; every section from 12 on needs ",
       "Vr_rot, ct_local and B_x.", call. = FALSE)
stopifnot(identical(as.integer(loc_12$best_pair), c(1L, 2L)))

pkg_pair <- as.integer(loc_op$best_pair)
Vr_B <- loc_12$Vr_rot                      # 4 x 2 -- the pinned (1,2) plane
Vr_C <- loc_op$Vr_rot                      # 4 x 2 -- the selected pair's plane
stopifnot(ncol(Vr_B) == 2L, ncol(Vr_C) == 2L,
          max(abs(crossprod(Vr_B) - diag(2))) < 1e-8,
          max(abs(crossprod(Vr_C) - diag(2))) < 1e-8)

cat(sprintf("Selected pair (bl_find_local_cf): (%d,%d)\n", pkg_pair[1], pkg_pair[2]))
cat("all_distances_mahalanobis (SQUARED -- the summary prints its square root):\n")
print(loc_op$all_distances_mahalanobis)

if (identical(pkg_pair, c(1L, 2L)))
  warning(sprintf(paste("Target row %d's selected pair IS (1,2), so views B and C are",
                        "the SAME search and the comparison shows nothing. Pick another",
                        "row -- 49 gives (2,3) on this data."), i), call. = FALSE)

prin_ang <- function(Q1, Q2)
  acos(pmin(1, pmax(-1, svd(t(Q1) %*% Q2)$d))) * 180 / pi

ang_BC <- prin_ang(Vr_B, Vr_C)
cat(sprintf("B vs C principal angles: %.2f, %.2f deg\n", ang_BC[1], ang_BC[2]))
if (ang_BC[2] < 5)
  warning(sprintf(paste("B and C planes are only %.2f deg apart -- the two views are",
                        "effectively the same picture. Pick a different target row."),
                  ang_BC[2]), call. = FALSE)

# Reproducibility check (note section 5.3). Cheap, and the alternative is quoting
# a counterfactual that a rerun does not reproduce.
reps <- replicate(5L, paste(as.integer(
  bl_find_local_cf(bl_res, tgt, verbose = FALSE)$best_pair), collapse = ","))
cat(sprintf("best_pair across 5 identical calls: %s\n",
            paste(unique(reps), collapse = " | ")))
if (length(unique(reps)) > 1L)
  warning(sprintf("best_pair is NOT reproducible for row %d: %s", i,
                  paste(unique(reps), collapse = " / ")), call. = FALSE)

# --- 7) .bl_rotate_full() — a mirror, for the per-pair table only -----------
# The package's .bl_rotate() (R/local_cf.R:45) is the authority for this maths.
# Two deliberate differences: this takes the FULL p x p V, and returns the FULL
# rotated basis, because section 9's table needs to know where every column landed
# whereas .bl_rotate() returns only the 2 a biplot needs.
#
# WHY A MIRROR AT ALL, given that section 6 already has the package's own planes.
# The section 9 table needs a rotation for EACH of the six candidate pairs, and
# bl_find_local_cf() returns Vr_rot only for the winner. There is no accessor for
# the losers' planes. So the mirror exists to fill that table and nothing else --
# NOTHING DRAWN IN THIS SCRIPT USES IT. Its faithfulness is asserted below against
# the two planes the package does expose.
#
# WHY THE PCs MUST BE bl_res$V AND NOT svd(Xstd)$v. Not a style choice. The
# Procrustes matrix is rank 1 (section 8), so the rotation is pinned on one
# direction and the rest is an arbitrary orthonormal completion -- and that
# completion is NOT stable to perturbation. bl_res$V and svd(Xstd)$v agree to
# ~4.65e-15, which sounds like nothing; rotating from the wrong one of them lands
# on a materially different plane. Since ct_local lives in Vr_rot's coordinates,
# drawing with a different plane would plot the package's contour into a frame it
# does not belong to. The measurement is printed below.
#
# WHAT IS AND IS NOT ASSERTED HERE -- this took a wrong turn on the first attempt
# and the distinction is the whole content of section 4.1.
#   ASSERTED: the mirror against the package's own private .bl_rotate() with
#     IDENTICAL inputs -- same x_target, same V, same pair. That is the only fair
#     test of whether the mirror implements the same maths, and it passes exactly.
#   REPORTED, NOT ASSERTED: the mirror against the package's LIVE objects
#     (loc_12$Vr_rot, loc_op$Vr_rot). bl_find_local_cf() reaches .bl_rotate()
#     through its own path, so its inputs differ from V_full at the 1e-15 level --
#     and section 8 measures that this is enough to land on a different plane. A
#     nonzero number in that column is the free completion, not an error. It comes
#     out LARGE for the pinned pair (1,2) and small for the selected pair, which is
#     exactly 17_1 section 15c's observation reproduced at p = 4.
# Nothing drawn in this script uses the mirror, so neither number affects a figure.
.bl_rotate_full <- function(x_target_st, V, proj_pair) {
  p    <- length(x_target_st)
  Vr   <- V[, proj_pair, drop = FALSE]                  # p x 2

  Y    <- rbind(-x_target_st, rep(0, p), x_target_st)   # 3 x p, RANK 1
  YV   <- Y %*% V                                       # 3 x p
  YVr  <- Y %*% Vr                                      # 3 x 2

  YVr_padded <- cbind(YVr, matrix(0, nrow = 3L, ncol = p - 2L))

  svd_res <- svd(t(YV) %*% YVr_padded)                  # p x p SVD
  A       <- svd_res$v %*% t(svd_res$u)                 # p x p rotation

  V %*% t(A)                                            # rotated full basis (p x p)
}

# Candidate pairs in the package's own order (.generate_pairs, colexicographic).
gen_pairs <- function(p) {
  out <- list()
  for (j in 2L:p) for (ii in 1L:(j - 1L)) out[[length(out) + 1L]] <- c(ii, j)
  do.call(rbind, out)
}
pairs_all  <- gen_pairs(p)
pair_names <- apply(pairs_all, 1L, function(r) sprintf("(%d,%d)", r[1], r[2]))

# Compare planes by PROJECTOR, never entry-wise: the in-plane basis is free to
# rotate, so entry-wise differences of order 1 are compatible with an identical
# plane. (CLAUDE.md section 3.)
Pj <- function(M) M %*% t(M)

cat("\n--- 7a) Mirror vs .bl_rotate(), IDENTICAL inputs — the fair test ---\n")
mirror_max <- 0
for (k in seq_len(nrow(pairs_all))) {
  pr   <- pairs_all[k, ]
  mine <- .bl_rotate_full(x_t, V_full, pr)[, 1:2, drop = FALSE]
  pkgr <- boundarylogic:::.bl_rotate(x_t, V_full, pr)$Vr_rot
  dd   <- max(abs(Pj(mine) - Pj(pkgr)))
  mirror_max <- max(mirror_max, dd)
  cat(sprintf("  pair %s plane diff: %.3e\n", pair_names[k], dd))
}
stopifnot(mirror_max < 1e-12)
cat("  -> .bl_rotate_full() reproduces .bl_rotate() exactly. The mirror is faithful.\n")

cat("\n--- 7b) Mirror vs the package's LIVE objects — reported, not asserted ---\n")
for (nm in c("pinned (1,2)", "selected")) {
  pr <- if (nm == "pinned (1,2)") c(1L, 2L) else pkg_pair
  vk <- if (nm == "pinned (1,2)") Vr_B else Vr_C
  m_pkg <- .bl_rotate_full(x_t, V_full,      pr)[, 1:2, drop = FALSE]
  m_svd <- .bl_rotate_full(x_t, svd(Xstd)$v, pr)[, 1:2, drop = FALSE]
  cat(sprintf(paste0("  %-13s rotated from bl_res$V : plane %.3e | entry-wise %.3e\n",
                     "  %-13s rotated from svd()$v  : plane %.3e\n"),
              nm, max(abs(Pj(m_pkg) - Pj(vk))), max(abs(m_pkg - vk)),
              "", max(abs(Pj(m_svd) - Pj(vk)))))
}
cat("  A nonzero plane difference here is the FREE COMPLETION, not an error: the\n")
cat("  package reaches .bl_rotate() through its own path, whose V differs from this\n")
cat("  script's at ~1e-15, and section 8b measures that this is enough. The gap is\n")
cat("  large for the pinned pair and small for the selected one -- 17_1 section 15c\n")
cat("  saw the same asymmetry at p = 3. Entry-wise differences are the in-plane spin.\n")
cat("  Nothing drawn uses the mirror; the drawn planes ARE loc_12/loc_op$Vr_rot.\n")

# --- 8) The rotation is under-determined — measured at p = 4 -----------------
# .bl_rotate() solves an orthogonal Procrustes problem against
# Y = rbind(-x, 0, x), whose rows are all multiples of x -- so Y is RANK 1, and so
# is the Procrustes matrix. The objective therefore imposes a single VECTOR
# constraint, A'u ~ w (u = V'x the target in PC coordinates, w = u with the
# non-pair entries zeroed), on a p x p rotation. The SVD of a rank-1 matrix has
# p - 1 exactly degenerate zero singular values, and the corresponding singular
# vectors -- hence A, hence the plane's second direction -- are whatever a
# particular LAPACK bidiagonalisation happens to emit for a degenerate subspace.
cat("\n--- 8) Under-determination ---\n")
u_t <- as.numeric(x_t %*% V_full)
for (k in seq_len(nrow(pairs_all))) {
  pr <- pairs_all[k, ]
  Yq <- rbind(-x_t, rep(0, p), x_t)
  Mq <- t(Yq %*% V_full) %*% cbind(Yq %*% V_full[, pr], matrix(0, 3L, p - 2L))
  d  <- svd(Mq)$d
  cat(sprintf("  pair %s Procrustes singular values: %s -> numerical rank %d of %d\n",
              pair_names[k], paste(format(d, digits = 3), collapse = "  "),
              sum(d > 1e-10 * max(d)), p))
  stopifnot(sum(d > 1e-10 * max(d)) == 1L)
}
cat("  Rank 1 of 4 -> a 3-dimensional null space. At p = 3 it was 1-dimensional.\n")

# ---- 8a) CAREFUL: free ROTATIONS and free PLANES are different counts -----
# The reference note's section 5.1 tabulates "subspace DOF ... minus contains the
# target ... free = 3" for p = 4, r = 2, and cites section 4.1's "3-dimensional
# manifold of rotations" as confirmation. The 3 is right FOR ROTATIONS and wrong
# FOR PLANES:
#   ROTATIONS. Those with A'u ~ w form a coset of Stab(u) = SO(p-1), so the
#     manifold has dim SO(4) - dim SO(3) = 6 - 3 = 3. This is section 4.1's number.
#   PLANES. 2-planes in 4-space containing a fixed line form Gr(1,3), of dimension
#     2 -- not 3. Containing a line is TWO conditions on Gr(2,4), not one.
# At p = 3 the note's "-1" is correct by accident (planes through a line in 3-space
# form a 1-parameter pencil), which is why 17_1 never tripped over it.
#
# The counts differ because a rotation fixes more than the plane: it also fixes
# where the target lands WITHIN the plane, leaving an in-plane spin -- exactly the
# entry-wise difference section 7 reports at a plane difference of ~1e-16.
#
# Measured rather than asserted. The tangent space to the family of 2-planes
# through the target is spanned by perturbing the plane's second direction inside
# the 2-dimensional orthogonal complement, so its dimension is read off the rank of
# the resulting derivative of the projector.
d0 <- Vr_B[, 2] - Vr_B[, 1] * sum(Vr_B[, 1] * Vr_B[, 2])
ray_u <- x_t / sqrt(sum(x_t^2))
d0 <- d0 - ray_u * sum(d0 * ray_u); d0 <- d0 / sqrt(sum(d0^2))
cmp <- svd(diag(p) - cbind(ray_u, d0) %*% t(cbind(ray_u, d0)))$u[, 1:(p - 2), drop = FALSE]
eps <- 1e-6
P0  <- Pj(cbind(ray_u, d0))
tng <- vapply(seq_len(ncol(cmp)), function(j) {
  dj <- d0 + eps * cmp[, j]; dj <- dj - ray_u * sum(dj * ray_u)
  dj <- dj / sqrt(sum(dj^2))
  as.numeric(Pj(cbind(ray_u, dj)) - P0) / eps
}, numeric(p * p))
sv_t <- svd(tng)$d
cat(sprintf("\n  8a) planes-through-the-target family: tangent singular values %s\n",
            paste(format(sv_t, digits = 3), collapse = "  ")))
cat(sprintf("      -> measured PLANE dimension = %d  (rotations = 3; the difference\n",
            sum(sv_t > 1e-8 * max(sv_t))))
cat("         is the in-plane spin, 3 = 2 + 1). The note's section 5.1 needs this fix.\n")
stopifnot(sum(sv_t > 1e-8 * max(sv_t)) == (p - 2L))

# ---- 8b) perturbation: section 4.1's table, on the live objects -----------
# Perturbing the input basis by amounts at or below ordinary floating-point noise
# produces materially different planes. The FIRST principal angle stays 0 -- the
# target's own direction is always shared -- so the SECOND is the arbitrary part.
#
# The re-orthonormalisation must stay NEAR V_full, or the table measures a sign
# flip instead of a perturbation. qr.Q() does not: it returns arbitrary column
# signs, so even at zero noise it lands ~1.8 away and every row reads the same.
# The polar factor U from the SVD is the NEAREST orthogonal matrix to M, which is
# the right tool and is what "re-orthonormalise only" means.
polar <- function(M) { s <- svd(M); s$u %*% t(s$v) }
cat("\n  8b) perturb the input basis, re-rotate onto the SAME pair, compare planes:\n")
set.seed(11L)
for (s in c(0, 1e-16, 1e-15, 1e-12, 1e-08)) {
  Vp <- if (s == 0) polar(V_full)
        else polar(V_full + matrix(rnorm(p * p, sd = s), p))
  vp <- .bl_rotate_full(x_t, Vp, c(1L, 2L))[, 1:2, drop = FALSE]
  v0 <- .bl_rotate_full(x_t, V_full, c(1L, 2L))[, 1:2, drop = FALSE]
  a  <- prin_ang(vp, v0)
  cat(sprintf("      perturbation %-8s | max|dV| %.2e | principal angles %6.2f, %6.2f deg\n",
              if (s == 0) "re-orth" else format(s), max(abs(Vp - V_full)), a[1], a[2]))
  # BOTH planes contain the target, so the first principal angle MUST be 0 -- and
  # at p = 4 that is real content, because two 2-planes in 4-space generically meet
  # only at the origin. This is the reference note's section 4.1 statement, and
  # this is the comparison it is about: perturbed versus unperturbed rotations of
  # ONE pair, not a rotated plane against the eigenplane its label names (section 9,
  # where the angle is emphatically NOT 0). The second angle is the arbitrary part.
  stopifnot(a[1] < 1e-6)
}

# ---- 8c) what that costs, WITHIN ONE PAIR --------------------------------
# The note's section 5.4 measures the consequence ACROSS pairs (pinned (1,2) versus
# the selected pair). The tighter claim -- that noise alone, at one fixed pair,
# moves the advice -- needs a within-pair measurement, and this is it.
#
# Two independent measurements, both of pair (1,2):
#   (a) the package itself, called twice. loc_12 (max_pairs = 1L) and the (1,2)
#       entry of loc_op's all_distances_mahalanobis are the SAME pair evaluated in
#       two calls. If the rotation were determined, they would be identical.
#   (b) the mirror, re-rotated from a basis perturbed at 1e-15, with the nearest
#       boundary point recomputed by the package's own recipe.
#
# plane_cf() reproduces bl_find_local_cf()'s within-plane search: a 200 x 200 grid
# over the padded range of the rotated training scores, mapped back to feature
# space, predicted, contoured at cutoff -/+ b_margin, wrong-side segments dropped,
# nearest surviving vertex by EUCLIDEAN distance in the plane, then that point
# mapped back and scored by SQUARED Mahalanobis. It omits the actionability filters
# (none are set) and the train_ranges test -- so it is a MEASUREMENT TOOL for this
# section only, not a fidelity check and not a drawing path.
plane_cf <- function(Vr, m = m_probe) {
  Zt  <- as.numeric(x_t %*% Vr)
  Ztr <- Xstd %*% Vr
  zr  <- range(Ztr); zpad <- diff(zr) * 0.1
  sq  <- seq(zr[1] - zpad, zr[2] + zpad, length.out = m)
  Zg  <- as.matrix(expand.grid(sq, sq))
  pg  <- pfun(Zg %*% t(Vr))
  ct  <- contourLines(sq, sq, matrix(pg, ncol = m, byrow = FALSE),
                      levels = c(cutoff - b_margin, cutoff + b_margin))
  if (!length(ct)) return(NULL)
  tgt_class <- as.integer(pfun(matrix(x_t, 1L)) >= cutoff)
  keep <- lapply(ct, function(s)
    if (as.integer(s$level >= cutoff) == tgt_class) NULL else cbind(s$x, s$y))
  keep <- keep[!vapply(keep, is.null, logical(1))]
  if (!length(keep)) return(NULL)
  Mi  <- do.call(rbind, keep)
  d2  <- (Mi[, 1] - Zt[1])^2 + (Mi[, 2] - Zt[2])^2
  idx <- which.min(d2)
  B_x <- as.numeric(Mi[idx, , drop = FALSE] %*% t(Vr))
  v   <- B_x - x_t
  list(B_x = B_x, dist_z = sqrt(d2[idx]),
       mah = if (is.null(metric_inv)) NA_real_
             else as.numeric(t(v) %*% metric_inv %*% v))
}

cat("\n  8c) the cost of the arbitrary completion, within pair (1,2):\n")
mah_12_pinned <- sqrt(loc_12$dist_mahalanobis)
mah_12_search <- sqrt(as.numeric(loc_op$all_distances_mahalanobis[["(1,2)"]]))
cat(sprintf(paste0("      (a) the package, called twice, same pair (1,2):\n",
                   "          max_pairs = 1L        sqrt(mah) = %.6f\n",
                   "          inside the full search sqrt(mah) = %.6f   difference %.2e\n"),
            mah_12_pinned, mah_12_search, abs(mah_12_pinned - mah_12_search)))
set.seed(23L)
Vq <- qr.Q(qr(V_full + matrix(rnorm(p * p, sd = 1e-15), p)))
cf0 <- plane_cf(.bl_rotate_full(x_t, V_full, c(1L, 2L))[, 1:2, drop = FALSE])
cf1 <- plane_cf(.bl_rotate_full(x_t, Vq,     c(1L, 2L))[, 1:2, drop = FALSE])
if (!is.null(cf0) && !is.null(cf1)) {
  cat(sprintf(paste0("      (b) the mirror, basis perturbed at 1e-15, same pair (1,2):\n",
                     "          B_x unperturbed (%s)\n",
                     "          B_x perturbed   (%s)\n",
                     "          the advice moves by %.4f in feature space; sqrt(mah) %.4f -> %.4f\n"),
              paste(sprintf("%+.3f", cf0$B_x), collapse = ", "),
              paste(sprintf("%+.3f", cf1$B_x), collapse = ", "),
              sqrt(sum((cf0$B_x - cf1$B_x)^2)), sqrt(cf0$mah), sqrt(cf1$mah)))
}
cat("      Reported, not asserted: the size of the shift is data-dependent, and\n")
cat("      the point is that it is nonzero at all from noise at 1e-15.\n")

# --- 9) What the pair label is worth ----------------------------------------
# Per candidate pair: theta (the rotation angle), the principal angles between the
# rotated plane and the eigenplane the pair NAMES, how much of the rotated plane's
# second direction actually lies inside that eigenplane, and which original PC each
# rotated column is closest to.
cat("\n--- 9) What the pair label is worth ---\n")
pl <- data.frame(pair = pair_names, theta_deg = NA_real_, ang1_deg = NA_real_,
                 ang2_deg = NA_real_, pct_in_named = NA_real_, row.names = NULL)
for (k in seq_len(nrow(pairs_all))) {
  pr <- pairs_all[k, ]
  Vd <- .bl_rotate_full(x_t, V_full, pr)
  ww <- u_t; ww[-pr] <- 0
  pl$theta_deg[k] <- acos(min(1, sqrt(sum(ww^2)) / sqrt(sum(u_t^2)))) * 180 / pi
  ang <- prin_ang(Vd[, 1:2, drop = FALSE], V_full[, pr, drop = FALSE])
  pl$ang1_deg[k] <- ang[1]; pl$ang2_deg[k] <- ang[2]
  tv <- as.numeric(x_t %*% Vd[, 1:2]); td <- Vd[, 1:2] %*% (tv / sqrt(sum(tv^2)))
  d2 <- Vd[, 1:2] %*% c(0, 1); d2 <- d2 - td * sum(d2 * td)
  d2 <- d2 / sqrt(sum(d2^2))
  Pn <- Pj(V_full[, pr, drop = FALSE])
  pl$pct_in_named[k] <- 100 * sum((Pn %*% d2)^2)
  cl <- function(col) which.min(acos(pmin(1, abs(as.numeric(t(col) %*% V_full)))))
  cat(sprintf("  pair %s: RC1 is closest to PC%d, RC2 to PC%d\n",
              pair_names[k], cl(Vd[, 1]), cl(Vd[, 2])))
}
print(pl, digits = 4)

# ---- 9a) A CORRECTION TO 17_1's SECTION 12 COMMENT -----------------------
# 17_1 section 12 says of this table: "the FIRST principal angle is 0 deg for every
# pair, but in 3D that is FORCED (any two 2-planes in 3-space share a line). The
# note reports the same 0 for p = 4, where it IS meaningful -- there two 2-planes
# generically meet only at the origin, so a 0 really does mark the target's
# direction."
#
# The first half is right; the parenthetical conflates two DIFFERENT comparisons,
# and this script is where that shows up.
#
#   THIS TABLE compares a rotated plane against the EIGENPLANE ITS LABEL NAMES.
#     The named eigenplane does not contain the target, so there is no reason for a
#     zero -- and measured at p = 4 there is none: ang1 runs from 1.11 to 55.62 deg
#     across the six pairs. At p = 3 it was 0 purely because two 2-planes in 3-space
#     must share a line.
#   THE NOTE'S SECTION 4.1 TABLE compares perturbed against unperturbed rotations
#     of ONE pair. BOTH of those contain the target, so ang1 = 0 necessarily, and at
#     p = 4 that is real content. That comparison is section 8b above, where the
#     zero is asserted, and section 12a, where the two drawn planes share the
#     target's ray to 12 decimal places.
#
# So the "0 is meaningful at p = 4" claim is true of 8b and 12a, and false of this
# table. ang1 here is a genuine measurement in its own right: how much of the
# target's direction the named eigenplane happens to capture.
cat(sprintf(paste0("\n  9a) ang1 here is NOT 0 at p = 4 (range %.2f to %.2f deg) -- see the\n",
                   "      comment above: this table compares against the NAMED eigenplane,\n",
                   "      which need not contain the target. The note's zero belongs to 8b.\n"),
            min(pl$ang1_deg), max(pl$ang1_deg)))
cat(sprintf(paste0("      Selected pair (%d,%d): its rotated plane's second direction lies\n",
                   "      %.1f%% inside the eigenplane its own label names.\n"),
            pkg_pair[1], pkg_pair[2],
            pl$pct_in_named[match(sprintf("(%d,%d)", pkg_pair[1], pkg_pair[2]), pl$pair)]))
cat("      (pct_in_named measures the second direction taken ORTHOGONAL to the\n")
cat("      target direction -- the definition is in the loop above, and a different\n")
cat("      convention gives a different number, so quote it with its definition.)\n")

# Does theta predict the overlap? The note's section 5.2 measures r = +0.33 at
# p = 3 and warns there is no cheap screen. Measured here at p = 4 rather than
# quoted, across every (row, pair) combination.
ov <- matrix(NA_real_, nrow(Xstd), nrow(pairs_all))
th <- ov
for (r0 in seq_len(nrow(Xstd))) {
  xr <- as.numeric(Xstd[r0, ]); ur <- as.numeric(xr %*% V_full)
  for (k in seq_len(nrow(pairs_all))) {
    pr <- pairs_all[k, ]
    Vd <- try(.bl_rotate_full(xr, V_full, pr), silent = TRUE)
    if (inherits(Vd, "try-error")) next
    wr <- ur; wr[-pr] <- 0
    th[r0, k] <- acos(min(1, sqrt(sum(wr^2)) / sqrt(sum(ur^2)))) * 180 / pi
    tv <- as.numeric(xr %*% Vd[, 1:2])
    if (sqrt(sum(tv^2)) < 1e-12) next
    td <- Vd[, 1:2] %*% (tv / sqrt(sum(tv^2)))
    d2 <- Vd[, 1:2] %*% c(0, 1); d2 <- d2 - td * sum(d2 * td)
    if (sqrt(sum(d2^2)) < 1e-12) next
    d2 <- d2 / sqrt(sum(d2^2))
    ov[r0, k] <- 100 * sum((Pj(V_full[, pr, drop = FALSE]) %*% d2)^2)
  }
}
ok <- is.finite(th) & is.finite(ov)
cat(sprintf(paste0("\n  theta vs %% inside the named eigenplane, over %d (row, pair)",
                   " combinations:\n    Pearson r = %+.2f | Spearman = %+.2f\n"),
            sum(ok), cor(th[ok], ov[ok]), cor(th[ok], ov[ok], method = "spearman")))
sm <- ok & th < 5
if (any(sm))
  cat(sprintf("    among the %d cases with theta < 5 deg the overlap spans %.1f%% to %.1f%%\n",
              sum(sm), min(ov[sm]), max(ov[sm])))
cat("  -> the overlap has to be computed; it cannot be inferred from theta.\n")

# --- 10) How the planes are tilted — principal angles, not one dihedral -----
# 17_1 section 17 writes a plane as s3 = s1*s1_ax + s2*s2_ax and reports the single
# dihedral tan^2(theta) = s1^2 + s2^2. That is the p = 3 case: the orthogonal
# complement is 1-dimensional, so there is exactly ONE nonzero principal angle.
#
# At p = 4 the complement is 2-dimensional, the tilt is a 2 x 2 MATRIX S mapping
# (PC1,PC2) coordinates to (PC3,PC4) coordinates, and there are TWO nonzero
# principal angles -- the arctangents of S's singular values. Derivation: with an
# orthonormal basis [I; S](I + S'S)^{-1/2} for the plane, the singular values of
# its inner product with [I; 0] are (1 + sigma_i^2)^{-1/2}, so cos(angle_i) =
# 1/sqrt(1 + sigma_i^2) and tan(angle_i) = sigma_i. Asserted numerically below
# against an independent principal-angle computation, not assumed.
#
# S versus the PC1-PC2 plane is generally FULL rank. Only the B-versus-C
# comparison has a forced zero first angle, because only those two planes share
# the target's ray.
cat("\n--- 10) Tilt of each drawn plane, relative to the PC1-PC2 plane ---\n")
tilt_report <- function(Vr, nm) {
  Q <- t(V_full) %*% Vr                        # 4 x 2, the plane in PC coordinates
  Q12 <- Q[1:2, , drop = FALSE]; Q34 <- Q[3:4, , drop = FALSE]
  a_ind <- prin_ang(Vr, V_full[, 1:2, drop = FALSE])
  if (abs(det(Q12)) < 1e-8) {
    cat(sprintf("  %-22s principal angles %6.2f, %6.2f deg  (S undefined: the plane is\n",
                nm, a_ind[1], a_ind[2]))
    cat("                         not transverse to span(PC3,PC4))\n")
    return(invisible(NULL))
  }
  S   <- Q34 %*% solve(Q12)
  sig <- svd(S)$d
  a_S <- atan(sig) * 180 / pi
  cat(sprintf(paste0("  %-22s principal angles %6.2f, %6.2f deg\n",
                     "                         atan(svd(S))    %6.2f, %6.2f deg  (match)\n",
                     "                         S = [%+.3f %+.3f ; %+.3f %+.3f]\n"),
              nm, a_ind[1], a_ind[2], sort(a_S)[1], sort(a_S)[2],
              S[1, 1], S[1, 2], S[2, 1], S[2, 2]))
  stopifnot(max(abs(sort(a_S) - sort(a_ind))) < 1e-6)
  invisible(S)
}
tilt_report(Vr_B, "B: pinned (1,2)")
tilt_report(Vr_C, sprintf("C: selected (%d,%d)", pkg_pair[1], pkg_pair[2]))
cat(sprintf("  B vs C:                principal angles %6.2f, %6.2f deg\n",
            ang_BC[1], ang_BC[2]))
cat("  The first B-vs-C angle is 0 because both planes contain the target (sec 12).\n")
cat("  A shallower tilt is easier to read; it is NOT a better answer -- selection is\n")
cat("  on the Mahalanobis distance to the counterfactual and nothing else.\n")

# --- 11) What the rotation costs -------------------------------------------
# 17_1 section 13's diagnostic block, at J = diag(1,1,0,0). Two differences from
# the 3-variable case:
#   - TWO components are dropped, not one, so the unrotated RSS equals S[3] + S[4];
#   - the unrotated basis is already lossy in 3D as well as 2D, which is the whole
#     subject of sections 16-17.
# J selects the 2D biplot: keep the first two axes, drop the rest.
sv_x   <- svd(Xstd)
D      <- sv_x$d; S_var <- D^2
stopifnot(min(D) > 1e-8)
J      <- diag(c(1, 1, 0, 0))
XtX    <- crossprod(Xstd)
tot_SS <- sum(diag(XtX))

# The drawn planes are 4 x 2; the diagnostics need a full p x p basis, so complete
# each with an orthonormal basis of its complement. adequacy indexes the ROWS of
# the basis, so it is PER VARIABLE, and J only ever selects columns 1-2 -- the
# completion affects the printed numbers only through those first two columns.
complete_basis <- function(Vr) {
  R  <- diag(p) - Pj(Vr)
  sv <- svd(R)
  cbind(Vr, sv$u[, 1:(p - 2), drop = FALSE])
}
V_A  <- V_full
V_Bf <- complete_basis(Vr_B)
V_Cf <- complete_basis(Vr_C)
stopifnot(max(abs(crossprod(V_Bf) - diag(p))) < 1e-8,
          max(abs(crossprod(V_Cf) - diag(p))) < 1e-8)

rc_diagnostics <- function(V_draw, label) {
  Uv  <- Xstd %*% V_draw %*% solve(diag(D))
  adq <- diag(V_draw %*% J %*% t(V_draw)) / diag(V_draw %*% t(V_draw))
  axp <- diag(V_draw %*% diag(S_var) %*% J %*% t(V_draw)) /
         diag(V_draw %*% diag(S_var) %*% t(V_draw))
  smp <- diag(Uv %*% diag(S_var) %*% J %*% t(Uv)) /
         diag(Uv %*% diag(S_var) %*% t(Uv))
  rss <- tot_SS - sum(diag(XtX %*% V_draw %*% J %*% t(V_draw)))
  xh  <- max(abs(x_t %*% V_draw[, 1:2] %*% t(V_draw[, 1:2]) - x_t))
  cat(sprintf(paste0("\n  %s\n",
      "    adequacy (per variable)  : %s\n",
      "    axis predictivity        : %s\n",
      "    RSS                      : %.4f  of total SS %.4f\n",
      "    target row %d sample pred : %.6f   | other %d points mean: %.4f\n",
      "    max|x_hat - x| (target)   : %.3e\n"),
      label, paste(sprintf("%.4f", adq), collapse = "  "),
      paste(sprintf("%.4f", axp), collapse = "  "),
      rss, tot_SS, i, smp[i], length(smp) - 1L, mean(smp[-i]), xh))
  list(adq = adq, axp = axp, smp = smp, rss = rss, xh = xh)
}

cat("\n--- 11) The cost of the rotation (J = diag(1,1,0,0): the 2D biplot) ---")
dg_A <- rc_diagnostics(V_A,  "A. Unrotated PCA")
dg_B <- rc_diagnostics(V_Bf, "B. RC biplot, rotated from (1,2)")
dg_C <- rc_diagnostics(V_Cf, sprintf("C. RC biplot, rotated from (%d,%d) [selected]",
                                     pkg_pair[1], pkg_pair[2]))

# adequacy is PER VARIABLE, so the cost is a REDISTRIBUTION across variables --
# some improve, some collapse -- not a uniform loss. Which variable collapses
# depends on where the rotated plane ends up, not on any ranking of the variables,
# so it is computed rather than assumed.
w_col <- which.min(dg_C$adq)
cat(sprintf(paste0("\n  THE TRADE, in one line (view C):\n",
    "    target row %d becomes EXACT       : sample predictivity %.4f -> %.6f\n",
    "    every other point gets worse      : mean %.4f -> %.4f\n",
    "    global residual rises             : RSS %.3f -> %.3f  (of %.0f)\n",
    "    representation is redistributed   : %s collapses, adequacy %.4f -> %.4f\n",
    "                                        (all variables: %s -> %s)\n"),
    i, dg_A$smp[i], dg_C$smp[i], mean(dg_A$smp[-i]), mean(dg_C$smp[-i]),
    dg_A$rss, dg_C$rss, tot_SS,
    feat_names[w_col], dg_A$adq[w_col], dg_C$adq[w_col],
    paste(sprintf("%.3f", dg_A$adq), collapse = "/"),
    paste(sprintf("%.3f", dg_C$adq), collapse = "/")))

stopifnot(
  abs(dg_A$rss - (S_var[3] + S_var[4])) < 1e-8,  # TWO components dropped, not one
  dg_B$rss > dg_A$rss,                           # PC1-PC2 minimises rank-2 residual
  dg_C$rss > dg_A$rss,
  abs(dg_B$smp[i] - 1) < 1e-8,                   # the target becomes exact ...
  abs(dg_C$smp[i] - 1) < 1e-8,
  dg_B$xh < 1e-8, dg_C$xh < 1e-8,                # ... x_hat = x for the target
  dg_A$xh > 1e-3                                 # ... and NOT before the rotation
)

# --- 12) W: a 3D basis that holds BOTH planes exactly -----------------------
# The construction, and the one-line proof that it exists. Plane C's residual after
# removing plane B has rank 1 -- C adds EXACTLY ONE direction to B -- so the two
# planes span 2 + 1 = 3 dimensions, and [Vr_B, that direction] is an orthonormal
# basis for their span.
#
# The frame is chosen deliberately: with W = [Vr_B, n3], plane B IS the z = 0
# plane, so the figure shows plane C tilting about the shared ray rather than two
# planes both at odd angles.
cat("\n--- 12) The exact 3D basis W = span(plane B, plane C) ---\n")
Rres <- Vr_C - Vr_B %*% (t(Vr_B) %*% Vr_C)
sv_R <- svd(Rres)
cat(sprintf("  residual of plane C after removing plane B: singular values %s\n",
            paste(format(sv_R$d, digits = 3), collapse = "  ")))
if (sv_R$d[1] <= 1e-8)
  stop("planes B and C coincide -- there is nothing to compare.", call. = FALSE)
if (sv_R$d[2] >= 1e-8)
  stop("planes B and C do not span exactly 3 dimensions; W would be invalid and ",
       "the 3D figure would be dishonest.", call. = FALSE)
cat("  -> rank 1: plane C adds EXACTLY ONE direction to plane B, so the two\n")
cat("     planes span exactly 3 dimensions and a 3D figure can be exact.\n")

n3 <- sv_R$u[, 1]
W  <- cbind(Vr_B, n3)                      # 4 x 3, orthonormal; plane B is z = 0

bx_B <- as.numeric(loc_12$B_x); bx_C <- as.numeric(loc_op$B_x)
in_W <- function(v) max(abs(as.numeric(v) - as.numeric(as.numeric(v) %*% W %*% t(W))))
stopifnot(
  max(abs(crossprod(W) - diag(3))) < 1e-10,          # orthonormal frame
  max(abs(Vr_B - W %*% t(W) %*% Vr_B)) < 1e-10,      # plane B lies in W
  max(abs(Vr_C - W %*% t(W) %*% Vr_C)) < 1e-10,      # plane C lies in W
  in_W(x_t) < 1e-10,                                 # the target lies in W
  in_W(bx_B) < 1e-10, in_W(bx_C) < 1e-10             # both counterfactuals do too
)
cat("  W is orthonormal; both planes, the target and both counterfactuals lie in\n")
cat("  it exactly. Only the 221 data points are projected (section 17).\n")

# ---- 12a) the shared ray, computed without assuming the answer -----------
# A shared direction is a nonzero (a; b) with Vr_B a = Vr_C b, i.e. a null vector
# of [Vr_B, -Vr_C]. At p = 4 a 1-dimensional null space is a RESULT: two 2-planes
# in 4-space generically meet only at the origin. At p = 3 it was forced.
Mint <- cbind(Vr_B, -Vr_C)
ns   <- svd(Mint)
cat(sprintf("\n  12a) [Vr_B, -Vr_C] singular values: %s -> %d zero\n",
            paste(format(ns$d, digits = 3), collapse = "  "), sum(ns$d < 1e-8)))
if (sum(ns$d < 1e-8) != 1L)
  stop("the two planes do not meet in exactly one line; the shared-ray argument ",
       "fails.", call. = FALSE)
v_int <- as.numeric(Vr_B %*% ns$v[1:2, p]); v_int <- v_int / sqrt(sum(v_int^2))
cos_ray <- abs(sum(v_int * ray_u))
cat(sprintf("  12a) the shared line IS the target's own ray: |cos| = %.12f\n", cos_ray))
stopifnot(abs(cos_ray - 1) < 1e-8)
cat("       Both planes contain the origin (a rotation has no translation) and the\n")
cat("       target, so their intersection is forced to be span(x_t).\n")

# --- 13) Slice diagnostics --------------------------------------------------
to_feature <- function(lo) {
  stopifnot(length(lo$ct_local) > 0L)
  M <- do.call(rbind, lapply(lo$ct_local, function(s) cbind(s$x, s$y)))
  list(xy = M, X = M %*% t(lo$Vr_rot))
}
s_B <- to_feature(loc_12)
s_C <- to_feature(loc_op)

# ---- 13a) both curves ARE p = 0.5 curves of the same 4D model -------------
cat("\n--- 13a) Both 2D boundaries are genuine p = 0.5 curves ---\n")
p_B <- pfun(s_B$X); p_C <- pfun(s_C$X)
cat(sprintf("  pinned (1,2)   : %4d pts | p in [%.5f, %.5f] | max|p-0.5| = %.2e\n",
            nrow(s_B$X), min(p_B), max(p_B), max(abs(p_B - 0.5))))
cat(sprintf("  selected (%d,%d) : %4d pts | p in [%.5f, %.5f] | max|p-0.5| = %.2e\n",
            pkg_pair[1], pkg_pair[2], nrow(s_C$X), min(p_C), max(p_C),
            max(abs(p_C - 0.5))))
stopifnot(max(abs(p_B - 0.5)) < 5e-3, max(abs(p_C - 0.5)) < 5e-3)
cat("  Neither is wrong -- they are different cross-sections of one 4D surface.\n")

# ---- 13b) ... but they are different curves ------------------------------
nn_dist <- function(A, B) apply(A, 1L, function(a) sqrt(min(colSums((t(B) - a)^2))))
d_bc <- nn_dist(s_B$X, s_C$X); d_cb <- nn_dist(s_C$X, s_B$X)
span <- max(dist(Xstd))
cat("\n--- 13b) They are genuinely DIFFERENT curves ---\n")
cat(sprintf(paste0("  nearest-neighbour distance: min %.4f | median %.4f | mean %.4f |\n",
                   "  max %.4f (Hausdorff %.4f), against a data span of %.3f\n"),
            min(d_bc), median(d_bc), mean(d_bc), max(d_bc),
            max(max(d_bc), max(d_cb)), span))
cat(sprintf("  -> they separate by up to %.0f%% of the data's extent.\n",
            100 * max(d_bc) / span))

# ---- 13c) the one point they must share ----------------------------------
cat("\n--- 13c) The one point they must share ---\n")
t_scan <- seq(0.02, 4, length.out = 40001)
p_ray  <- pfun(outer(t_scan, x_t))
k0     <- which.min(abs(p_ray - 0.5))
crossed <- abs(p_ray[k0] - 0.5) < 0.01
n_cross <- sum(abs(diff(sign(p_ray - 0.5))) > 0)
if (!crossed) {
  cat("  the surface does NOT meet the ray within the scanned range -- no shared point.\n")
  x_cross <- NULL
} else {
  if (k0 == 1L || k0 == length(t_scan))
    warning("ray crossing sits at the scan boundary; widen the t range.", call. = FALSE)
  x_cross <- t_scan[k0] * x_t
  # k0 is the scan point CLOSEST to 0.5, which need not be the FIRST crossing --
  # the surface can meet the ray more than once. n_cross says how many times it
  # does for the target actually set.
  cat(sprintf(paste0("  surface crosses the ray %d time(s); closest scan point at\n",
                     "  t = %.5f -> (%s), p = %.6f\n"),
              n_cross, t_scan[k0],
              paste(sprintf("%+.4f", x_cross), collapse = ", "), p_ray[k0]))
  cat(sprintf("  (the target is at t = 1, so the crossing lies %s it)\n",
              if (t_scan[k0] < 1) "between the origin and" else "beyond"))
  dB <- min(sqrt(colSums((t(s_B$X) - x_cross)^2)))
  dC <- min(sqrt(colSums((t(s_C$X) - x_cross)^2)))
  step <- diff(loc_12$xseq)[1]
  cat(sprintf("  distance from that point to each curve: pinned %.5f | selected %.5f\n",
              dB, dC))
  cat(sprintf("  contour grid step, for scale: %.5f -> both pass through it\n", step))
  stopifnot(dB < 3 * step, dC < 3 * step)
}

# ---- 13d) the in-plane axes differ too (the cosmetic second reason) ------
cat("\n--- 13d) The panels are also drawn on different in-plane axes ---\n")
cat(sprintf("  RC1 differs by %.2f deg | RC2 by %.2f deg\n",
            acos(min(1, abs(sum(Vr_B[, 1] * Vr_C[, 1])))) * 180 / pi,
            acos(min(1, abs(sum(Vr_B[, 2] * Vr_C[, 2])))) * 180 / pi))
cat(sprintf("  target at (%.4f, %.4f) vs (%.4f, %.4f) -- norms %.4f and %.4f, both = ||x_t||\n",
            loc_12$Z_target[1], loc_12$Z_target[2],
            loc_op$Z_target[1], loc_op$Z_target[2],
            sqrt(sum(loc_12$Z_target^2)), sqrt(sum(loc_op$Z_target^2))))
stopifnot(abs(sqrt(sum(loc_12$Z_target^2)) - sqrt(sum(x_t^2))) < 1e-8,
          abs(sqrt(sum(loc_op$Z_target^2)) - sqrt(sum(x_t^2))) < 1e-8)
cat("  So even shared structure appears rotated between the two panels. THIS part\n")
cat("  is cosmetic -- the norms confirm both planes contain the target.\n")

# ---- 13e) why it is NOT cosmetic ----------------------------------------
cat("\n--- 13e) Why it is NOT cosmetic ---\n")
for (nm in c("pinned (1,2)", "selected")) {
  lo <- if (nm == "pinned (1,2)") loc_12 else loc_op
  cat(sprintf("  %-13s B_x = (%s) | dist_z %.4f | sqrt(mah) %.4f\n", nm,
              paste(sprintf("%+.3f", as.numeric(lo$B_x)), collapse = ", "),
              lo$dist_z, sqrt(lo$dist_mahalanobis)))
}
cat(sprintf("  the two counterfactuals are %.4f apart in feature space.\n",
            sqrt(sum((bx_B - bx_C)^2))))
cat("  Different slice -> different nearest boundary point -> different advice.\n")

# --- 14) The 2D biplots — the boundary drawn only where it is exact ---------
# One code path, three calls, so the panels differ in NOTHING but the basis and
# the contour. Axis labels: PC1/PC2 for A, which is genuinely unrotated; RC1/RC2
# for B and C. The named pair is provenance for WHICH PLANE was rotated, so it goes
# in the title, never on an axis.
make_2d_biplot <- function(Vr2, ct, ttl, axis_lab, cf_col = NULL, cf = NULL) {
  Z <- as.data.frame(Xstd %*% Vr2); names(Z) <- c("s1", "s2")
  z_t <- as.numeric(x_t %*% Vr2)

  f <- plot_ly() |>
    add_markers(data = Z, x = ~s1, y = ~s2,
                marker = list(size = 7, color = pt_col, opacity = 0.4),
                name = "Data points", showlegend = FALSE)
  for (k in seq_len(nrow(Vr2)))
    f <- f |> add_trace(x = c(0, Vr2[k, 1]) * 3, y = c(0, Vr2[k, 2]) * 3,
                        type = "scatter", mode = "lines+text",
                        text = c("", feat_names[k]), textposition = "top center",
                        textfont = list(size = 14, color = "#222222"),
                        line = list(color = "#5E88FC", width = 5), showlegend = FALSE)
  for (s in ct)
    f <- f |> add_trace(x = s$x, y = s$y, type = "scatter", mode = "lines",
                        line = list(width = 6, color = "black"),
                        name = "Decision boundary", showlegend = FALSE)
  if (!is.null(cf)) {
    z_cf <- as.numeric(cf %*% Vr2)
    f <- f |>
      add_trace(x = c(z_t[1], z_cf[1]), y = c(z_t[2], z_cf[2]),
                type = "scatter", mode = "lines",
                line = list(color = cf_col, width = 3, dash = "dot"),
                showlegend = FALSE) |>
      add_markers(x = z_cf[1], y = z_cf[2],
                  marker = list(size = 11, symbol = "circle", color = cf_col),
                  name = "Counterfactual", showlegend = FALSE)
  }
  f |>
    add_markers(x = z_t[1], y = z_t[2],
                marker = list(size = 13, symbol = "diamond", color = "magenta"),
                name = "Target", showlegend = FALSE) |>
    layout(title = ttl,
           xaxis = list(title = axis_lab[1], scaleanchor = "y", scaleratio = 1),
           yaxis = list(title = axis_lab[2]), showlegend = FALSE)
}

# View A's boundary is the GLOBAL biplot's, in PC1-PC2 with PC3 and PC4 pinned at
# 0. It is the population's boundary, not this target's -- the target does not lie
# in that plane, which is the whole motivation for rotating.
bip_A <- make_2d_biplot(
  V_full[, 1:2], bl_res$biplot_grid$ct,
  "A. Unrotated PCA biplot - the POPULATION's boundary (the target is off-plane)",
  c("PC1", "PC2"))
bip_B <- make_2d_biplot(
  Vr_B, loc_12$ct_local,
  "B. RC biplot - rotated from the (1,2) plane [pinned]",
  c("RC1", "RC2"), "red", bx_B)
bip_C <- make_2d_biplot(
  Vr_C, loc_op$ct_local,
  sprintf("C. RC biplot - rotated from the (%d,%d) plane [selected]",
          pkg_pair[1], pkg_pair[2]),
  c("RC1", "RC2"), "darkblue", bx_C)

# The rotation's whole job: put the target ON the drawn plane. Verified on the
# residual outside each plane -- exactly 0 for B and C, nonzero for A.
off_A <- sqrt(sum((x_t - as.numeric(x_t %*% V_full[, 1:2] %*% t(V_full[, 1:2])))^2))
off_B <- sqrt(sum((x_t - as.numeric(x_t %*% Vr_B %*% t(Vr_B)))^2))
off_C <- sqrt(sum((x_t - as.numeric(x_t %*% Vr_C %*% t(Vr_C)))^2))
cat(sprintf("\n--- 14) Target's distance from each drawn plane ---\n"))
cat(sprintf("  A (PC1-PC2): %.4f   B (pinned): %.3e   C (selected): %.3e\n",
            off_A, off_B, off_C))
stopifnot(off_A > 1e-6, off_B < 1e-8, off_C < 1e-8)

p_glob <- pfun(do.call(rbind, lapply(bl_res$biplot_grid$ct,
                                     function(s) cbind(s$x, s$y))) %*% t(V_full[, 1:2]))
cat(sprintf("  view A's global boundary: max|p - 0.5| = %.2e (it too is exact)\n",
            max(abs(p_glob - 0.5))))
stopifnot(max(abs(p_glob - 0.5)) < 5e-3)

# --- 15/16) One 3D code path, two bases -------------------------------------
# make_3d() draws the same objects in whatever 3-column basis it is given, so the
# exact figure (section 15, basis W) and the failure exhibit (section 16, basis
# PC1-PC2-PC3) differ in NOTHING but the basis. That is the only way the comparison
# is trustworthy.
#
# NO DECISION SURFACE IS DRAWN. The p = 0.5 locus of a 4-variable model is a
# 3-dimensional hypersurface in standardised feature space; no 2D surface in a 3D
# box is that object, and the previous version's smoothed green sheet was a
# regression through a probability band, not a level set. What is drawn instead is
# exact: the two planes, and on each the package's own boundary contour.
make_3d <- function(B3, ttl, exact) {
  Zb   <- Xstd %*% B3
  rad  <- max(sqrt(rowSums(Zb^2)))
  toB  <- function(X) X %*% B3                       # feature space -> this basis
  a_t  <- as.numeric(x_t %*% B3)

  quad <- function(Vr) {
    e1 <- as.numeric(t(B3) %*% Vr[, 1]); e2 <- as.numeric(t(B3) %*% Vr[, 2])
    L1 <- max(max(abs(Zb %*% e1)), 0.15 * rad) * 1.05
    L2 <- max(max(abs(Zb %*% e2)), 0.15 * rad) * 1.05
    rbind(-L1 * e1 - L2 * e2, L1 * e1 - L2 * e2,
           L1 * e1 + L2 * e2, -L1 * e1 + L2 * e2)
  }
  add_quad <- function(fig, cor, col, nm)
    add_trace(fig, type = "mesh3d", x = cor[, 1], y = cor[, 2], z = cor[, 3],
              i = c(0, 0), j = c(1, 2), k = c(2, 3),
              facecolor = rep(col, 2), opacity = 0.22,
              hoverinfo = "skip", name = nm, showlegend = TRUE)

  # The counterfactual legend carries the MAHALANOBIS distance, because that -- and
  # nothing else -- is what "selected" means. dist_mahalanobis is stored SQUARED,
  # so its square root is the distance itself.
  cf_lab <- function(lo, nm) {
    if (is.null(lo$B_x)) return(NA_character_)
    if (is.na(lo$dist_mahalanobis))
      sprintf("CF | %s - Euclidean %.4f (Mahalanobis unavailable)", nm, lo$dist_z)
    else sprintf("CF | %s - Mahalanobis %.4f", nm, sqrt(lo$dist_mahalanobis))
  }
  sel_B <- !is.na(loc_12$dist_mahalanobis) && !is.na(loc_op$dist_mahalanobis) &&
           loc_12$dist_mahalanobis < loc_op$dist_mahalanobis
  lab_B <- paste0(cf_lab(loc_12, "pinned (1,2)"),
                  if (isTRUE(sel_B)) "  <- selected (smallest)" else "")
  lab_C <- paste0(cf_lab(loc_op, sprintf("selected (%d,%d)", pkg_pair[1], pkg_pair[2])),
                  if (!isTRUE(sel_B)) "  <- selected (smallest)" else "")

  add_cf <- function(fig, bx, col, lab) {
    b <- as.numeric(bx %*% B3)
    fig |>
      add_trace(x = c(a_t[1], b[1]), y = c(a_t[2], b[2]), z = c(a_t[3], b[3]),
                type = "scatter3d", mode = "lines",
                line = list(color = col, width = 4, dash = "dot"),
                hoverinfo = "skip", showlegend = FALSE) |>
      add_markers(x = b[1], y = b[2], z = b[3],
                  marker = list(size = 9, color = col, symbol = "circle"),
                  name = lab, showlegend = TRUE)
  }

  cB <- toB(s_B$X); cC <- toB(s_C$X)
  ray <- as.numeric(ray_u %*% B3)
  loss <- pmax(0, 1 - rowSums(Zb^2) / rowSums(Xstd^2))

  fig <- plot_ly() |>
    add_markers(x = Zb[, 1], y = Zb[, 2], z = Zb[, 3],
                marker = list(size = 5, color = pt_col, opacity = 0.4),
                name = sprintf("Data points (projected; median loss %.3f)",
                               median(loss)),
                showlegend = TRUE)
  for (k in seq_len(nrow(B3)))
    fig <- fig |> add_trace(x = c(0, B3[k, 1]) * 3, y = c(0, B3[k, 2]) * 3,
                            z = c(0, B3[k, 3]) * 3,
                            type = "scatter3d", mode = "lines+text",
                            text = c("", feat_names[k]), textposition = "top center",
                            textfont = list(size = 14, color = "#222222"),
                            line = list(color = "#5E88FC", width = 5),
                            showlegend = FALSE)
  fig <- fig |>
    add_quad(quad(Vr_B), "purple", "plane: pinned (1,2)") |>
    add_quad(quad(Vr_C), "blue",
             sprintf("plane: selected (%d,%d)", pkg_pair[1], pkg_pair[2])) |>
    add_trace(x = c(-rad, rad) * ray[1], y = c(-rad, rad) * ray[2],
              z = c(-rad, rad) * ray[3], type = "scatter3d", mode = "lines",
              line = list(color = "#666666", width = 4, dash = "dash"),
              name = "target's ray (where the planes meet)", showlegend = TRUE) |>
    add_trace(x = cB[, 1], y = cB[, 2], z = cB[, 3],
              type = "scatter3d", mode = "lines",
              line = list(color = "red", width = 6),
              name = "boundary: pinned (1,2)", showlegend = TRUE) |>
    add_trace(x = cC[, 1], y = cC[, 2], z = cC[, 3],
              type = "scatter3d", mode = "lines",
              line = list(color = "darkblue", width = 6),
              name = sprintf("boundary: selected (%d,%d)", pkg_pair[1], pkg_pair[2]),
              showlegend = TRUE)
  fig <- add_cf(fig, bx_B, "red",      lab_B)
  fig <- add_cf(fig, bx_C, "darkblue", lab_C)
  if (!is.null(x_cross)) {
    xc <- as.numeric(x_cross %*% B3)
    fig <- fig |>
      add_markers(x = xc[1], y = xc[2], z = xc[3],
                  marker = list(size = 10, color = "black", symbol = "circle"),
                  name = if (exact) "both curves pass here"
                         else "where both curves SHOULD pass", showlegend = TRUE)
  }
  fig |>
    add_markers(x = a_t[1], y = a_t[2], z = a_t[3],
                marker = list(size = 9, symbol = "diamond", color = "magenta"),
                name = "Target", showlegend = TRUE) |>
    layout(title = ttl,
           scene = list(aspectmode = "data",
                        xaxis = list(title = if (exact) "W1 (= RC1 of plane B)" else "PC1"),
                        yaxis = list(title = if (exact) "W2 (= RC2 of plane B)" else "PC2"),
                        zaxis = list(title = if (exact) "W3 (the direction C adds)" else "PC3"),
                        camera = list(eye = list(x = 1.4, y = 1.4, z = 0.9))),
           showlegend = TRUE)
}

fig_exact <- make_3d(
  W, paste("Both rotated planes on ONE exact 3D basis - each 2D biplot is a",
           "cross-section, not a view.<br>",
           "<sup>Planes, boundaries, target and both counterfactuals are EXACT here;",
           "only the data points are projected. Selection is by smallest",
           "Mahalanobis distance.</sup>"), exact = TRUE)

fig_broken <- make_3d(
  V_full[, 1:3],
  paste("THE SAME FIGURE IN THE PCA BASIS - METRICALLY UNTRUSTWORTHY.<br>",
        "<sup>Dropping PC4 preserves what touches what, but foreshortens every",
        "distance and angle - the two boundaries look more alike, and the",
        "counterfactuals closer, than they are. See section 16.</sup>"),
  exact = FALSE)

# --- 16) Where the 3D picture of 4D data fails — measured -------------------
# BE PRECISE ABOUT THE FAILURE MODE. An orthogonal projection is LINEAR, so it
# preserves incidence: a 2-plane through the origin projects to a 2-plane through
# the origin, a curve lying in a plane projects to a curve lying in the projected
# plane, and the point the two curves share projects to a point the two projected
# curves still share. A first attempt at this section measured the shared point's
# distance to each curve before and after projecting and found it essentially
# unchanged (0.0166 -> 0.0160) -- because it MUST be. Nothing about the picture's
# topology breaks, and claiming it does would be wrong.
#
# What projection destroys is the METRIC. Every distance and angle is foreshortened
# by an amount that depends on direction, so the reader cannot trust anything they
# measure off the figure: how far the counterfactual is from the target, how far
# apart the two boundaries are, or how steeply the planes cross. Foreshortening is
# also one-sided -- it can only shrink -- so two curves that are far apart in 4D
# can be made to look nearly coincident, and never the reverse.
#
# In W none of this arises: every object the argument depends on lies IN W, and an
# orthonormal basis restricted to its own span is an isometry, so all of those
# distances and angles are exact.
cat("\n--- 16) Where a 3D picture of 4D data fails: the metric, not the topology ---\n")
pc4   <- V_full[, 4]
projA <- V_full[, 1:3]

cat(sprintf(paste0("  what the PCA view discards, as a PC4 component:\n",
                   "    pinned (1,2) curve   max|. PC4| = %.4f\n",
                   "    selected (%d,%d) curve max|. PC4| = %.4f\n",
                   "    the target           |. PC4| = %.4f  (and PC3 = %+.4f)\n",
                   "    the two planes       max|direction . PC4| = %.4f / %.4f\n",
                   "  against a data span of %.3f -- so up to %.0f%% of the extent is\n",
                   "  projected away from objects the argument depends on.\n"),
            max(abs(s_B$X %*% pc4)), pkg_pair[1], pkg_pair[2],
            max(abs(s_C$X %*% pc4)), abs(sum(x_t * pc4)), sum(x_t * V_full[, 3]),
            max(abs(t(Vr_B) %*% pc4)), max(abs(t(Vr_C) %*% pc4)),
            span, 100 * max(abs(s_C$X %*% pc4)) / span))

# Foreshortening, measured on the quantities a reader would actually read off.
d3 <- function(a, b, B3) sqrt(sum((as.numeric(a) %*% B3 - as.numeric(b) %*% B3)^2))
d4 <- function(a, b) sqrt(sum((as.numeric(a) - as.numeric(b))^2))
cat("\n  16a) foreshortening of the quantities a reader measures off the figure:\n")
cat(sprintf("  %-34s %9s %9s %9s %9s\n", "", "true", "in W", "in PCA", "PCA/true"))
rows <- list(
  list("target -> CF, pinned (1,2)",  x_t, bx_B),
  list("target -> CF, selected",      x_t, bx_C),
  list("CF to CF",                    bx_B, bx_C))
for (r0 in rows) {
  tv <- d4(r0[[2]], r0[[3]])
  cat(sprintf("  %-34s %9.4f %9.4f %9.4f %9.3f\n", r0[[1]], tv,
              d3(r0[[2]], r0[[3]], W), d3(r0[[2]], r0[[3]], projA),
              d3(r0[[2]], r0[[3]], projA) / tv))
}
# curve separation: the headline number of section 13b
nn_in <- function(B3) {
  A <- s_B$X %*% B3; Bm <- s_C$X %*% B3
  max(apply(A, 1L, function(a) sqrt(min(colSums((t(Bm) - a)^2)))))
}
cat(sprintf("  %-34s %9.4f %9.4f %9.4f %9.3f\n", "max curve-to-curve separation",
            max(d_bc), nn_in(W), nn_in(projA), nn_in(projA) / max(d_bc)))
# the angle between the two planes, as drawn
ang_in <- function(B3) {
  qB <- qr.Q(qr(t(B3) %*% Vr_B)); qC <- qr.Q(qr(t(B3) %*% Vr_C))
  prin_ang(qB, qC)[2]
}
cat(sprintf("  %-34s %8.2f%s %8.2f%s %8.2f%s %9.3f\n", "plane B vs C, 2nd angle (deg)",
            ang_BC[2], "", ang_in(W), "", ang_in(projA), "",
            ang_in(projA) / ang_BC[2]))
cat("\n  W reproduces every row exactly (ratio 1.000): the objects lie in it, and an\n")
cat("  orthonormal basis is an isometry on its own span.\n")
cat("  DISTANCES can only shrink under projection, and the worst offender is the one\n")
cat("  that matters most -- the target-to-counterfactual distance, i.e. how much\n")
cat("  would have to change for the decision to flip. It reads short by up to\n")
cat(sprintf("  %.0f%% in the PCA figure.\n",
            100 * (1 - min(d3(x_t, bx_B, projA) / d4(x_t, bx_B),
                           d3(x_t, bx_C, projA) / d4(x_t, bx_C)))))
cat("  ANGLES are not distances and can go EITHER way: the plane-to-plane angle\n")
cat("  comes out slightly LARGER as drawn than it is. So 'the projection shrinks\n")
cat("  everything' would be wrong -- what it does is distort, direction by\n")
cat("  direction, and nothing measured off the figure can be trusted either way.\n")
cat("  What it does NOT do is break incidence: the shared point stays shared. The\n")
cat("  figure is metrically untrustworthy, not topologically broken.\n")

# --- 17) The loss budget, both bases side by side ---------------------------
# The ranking inverts, and that is the point. PCA retains more of the DATA CLOUD --
# that is exactly what it is built to do -- and destroys the GEOMETRY. Neither
# basis is "better" without naming what it is better at.
cat("\n--- 17) Loss budget: W versus PC1-PC2-PC3 ---\n")
loss_of <- function(B3) pmax(0, 1 - rowSums((Xstd %*% B3)^2) / rowSums(Xstd^2))
lW <- loss_of(W); lA <- loss_of(projA)
cat(sprintf("  %-32s %16s %16s\n", "", "W", "PC1-PC2-PC3"))
cat(sprintf("  %-32s %16s %16s\n", "planes stay planes",
            "yes (contained)", "yes (projected)"))
cat(sprintf("  %-32s %16s %16s\n", "curves stay in their planes", "yes", "yes"))
cat(sprintf("  %-32s %16s %16s\n", "shared point stays shared", "yes", "yes"))
cat(sprintf("  %-32s %16s %16s\n", "distances (target -> CF)",
            "exact", sprintf("to %.0f%% short",
              100 * (1 - min(d3(x_t, bx_B, projA) / d4(x_t, bx_B),
                             d3(x_t, bx_C, projA) / d4(x_t, bx_C))))))
cat(sprintf("  %-32s %16s %16s\n", "angle between the planes",
            "exact", sprintf("%.2f vs %.2f deg", ang_in(projA), ang_BC[2])))
cat(sprintf("  %-32s %16s %16s\n", "geometry's discarded component",
            "none", sprintf("up to %.3f", max(abs(s_C$X %*% pc4)))))
cat(sprintf("  %-32s %16s %16s\n", "data points: median lost frac",
            sprintf("%.4f", median(lW)), sprintf("%.4f", median(lA))))
cat(sprintf("  %-32s %16s %16s\n", "data points: max lost frac",
            sprintf("%.4f", max(lW)), sprintf("%.4f", max(lA))))
cat("  -> W is exact for every object the argument depends on, and keeps the data\n")
cat("     cloud worse. PCA keeps the cloud better, and foreshortens every\n")
cat("     measurement the reader would take off the geometry. Choose the basis by\n")
cat("     what the figure is FOR: here, comparing two slices, so W.\n")

# --- 18) Figures ------------------------------------------------------------
# Run these one at a time in RStudio.

# The three 2D biplots -- the only place a decision boundary is drawn, and the
# only place it is exact without qualification.
bip_A
bip_B
bip_C

# THE HEADLINE FIGURE. Two translucent planes cutting through one cloud, each
# carrying its own boundary curve, the dashed grey ray where the planes meet, and
# the black marker where both curves pass through it. Everything geometric is
# exact. Rotate it: the two curves are visibly different objects, and they touch
# only on the ray.
fig_exact

# THE FAILURE EXHIBIT. The identical construction in the PCA basis. It still LOOKS
# coherent -- and that is the trap. Incidence survives projection, so nothing looks
# obviously broken; what is wrong is every distance and angle in it, foreshortened
# by up to the ratios section 16 tabulates. Shown to be measured, not trusted.
fig_broken

# The package's own biplots, for side-by-side comparison. plot(loc_12) should match
# bip_B's boundary and plot(loc_op) should match bip_C's, since both take ct_local
# verbatim; plot(bl_res) should match bip_A's.
plot(bl_res)
plot(loc_12)
plot(loc_op)
