# Originated in the `Explain XAI code` project's sandbox as
# `R/In-development/17_1_loan-gam-pca-3d-3var-rotation.R` (snapshot 2026-09-10).
# This repo is now this script's home; edit it here.

# ============================================================================
# The rotation, demonstrated losslessly — GAM decision surface, loan data, 3 vars
# ----------------------------------------------------------------------------
# WHY THIS SCRIPT EXISTS
#
# Script 17_2 shows the same rotation on FOUR features reduced to 3 PCs. There,
# two separate effects are tangled together:
#   1. the rotation itself  — a change of frame; nothing is lost, and
#   2. the 4D -> 3D truncation — genuinely lossy for every point but the target.
# A reader cannot tell which effect produces what they see.
#
# This script drops to THREE features, so PCA gives a 3 x 3 orthogonal V and the
# 3D view is EXACTLY LOSSLESS FOR EVERY POINT (asserted in section 4). Anything
# that differs between the three views below is therefore attributable to the
# rotation alone. That is the whole point: it makes the rotation confirmable.
#
# THE THREE VIEWS (all built by one function, so only the basis differs)
#   A. Unrotated        — the target FLOATS OFF the slot-3 = 0 plane, so the 2D
#                         biplot's boundary is not the target's boundary. Motivation.
#   B. Rotated -> (1,2) — the target lands exactly ON the plane.
#   C. Rotated -> optimal pair — the plane the package actually uses (here (1,3)).
#
# Each view is drawn twice: as a 3D surface and as the flat 2D biplot you get by
# looking down the third axis. The black curve is the surface's slot-3 = 0
# contour — the decision boundary the 2D biplot draws.
#
# Section 9 replicates the package's pair-selection rule and prints it beside the
# package's own numbers; section 12 checks the rotation against boundarylogic.
# Mirrors script 15 (same data, same GAM) plus the rotation.
#
# ONE RESULT WORTH KNOWING UP FRONT (section 12c). The rotated plane is NOT
# uniquely determined. .bl_rotate() solves an orthogonal Procrustes problem whose
# matrix is RANK 1, because Y = rbind(-x, 0, x) has rank 1. So the rotation is
# pinned on one direction only; the rest is an arbitrary orthonormal completion,
# and ANY plane containing the target is a valid answer. What IS guaranteed is
# that the target lies in the plane and keeps its in-plane coordinates.
# The practical consequence: bl_res$V differs from svd(Xstd)$v at ~1e-15, and
# that is enough to land on a different plane -- it does for pair (1,2) and does
# not for (1,3). Since the nearest boundary point depends on the plane, dist_z
# and the squared Mahalanobis distance inherit the sensitivity. So section 12
# asserts the mirror against .bl_rotate() with IDENTICAL inputs (exact, diff 0),
# and merely REPORTS the comparison against the package's live objects.
#
# ----------------------------------------------------------------------------
# HOW TO READ THE ROTATED VIEWS: THEY ARE **RC BIPLOTS**, NOT "ROTATED PCs"
#
# After the rotation the axes are RC1 / RC2 -- "rotated components", the thesis's
# own term for a rotated component (RC) biplot. They are NOT "rotated PC1 / PC3",
# and the named pair is NOT an axis identity. From the reference note
# (biplot_rotation_and_pair_selection.md, section 2.7):
#
#   "the SVD construction guarantees the target's rotated information concentrates
#    in columns 1-2 of Vrho REGARDLESS of which eigenvector pair was used to build
#    Vr -- which is why the function can simply always return columns 1-2 rather
#    than tracking where the pair's original indices ended up."
#
# SECTION 12 PRINTS THE EVIDENCE for whatever target is set -- per pair: the
# rotation angle theta, the principal angles against the eigenplane the pair NAMES,
# how much of the rotated plane actually lies inside that eigenplane, and which
# original PC each rotated column is closest to. Read it there; no numbers are
# hardcoded in this header, so it does not go stale when the target changes.
#
# Two things that table reliably shows:
#   - the pair label does not survive into the axes. It is common for the rotated
#     plane to lie mostly OUTSIDE the eigenplane its label names, and for the
#     columns not even to preserve the pair's order;
#   - the FIRST principal angle is 0 deg for every pair, but in 3D that is FORCED:
#     any two 2-planes in 3-space intersect in a line. (The note reports the same 0
#     for p = 4, where it IS meaningful -- there two 2-planes generically meet only
#     at the origin, so a 0 really does mark the target's direction. It does not
#     carry over here.) The SECOND angle is the informative one.
# The pair selects WHICH PLANE gets rotated; it is not an axis identity.
#
# THE LINK TO THE 1D EXAMPLE ("A3 Rotated Biplot.R")
#
# A3 reduces 2 variables to a ONE-dimensional biplot (r = 1). After rotating onto
# PC1 only that single component remains, and PC2 simply disappears -- there is no
# "rotated PC2" to speak of. This script is the same construction at r = 2: after
# rotation only RC1 and RC2 remain, whatever pair was fed in, and the third axis is
# a residual the 2D biplot drops. Its adequacy collapsing from 0.9999 to 0.072
# (section 13) is the measurable trace of that disappearance.
# ============================================================================

rm(list = ls())

library(boundarylogic)
if (!requireNamespace("plotly",   quietly = TRUE)) install.packages("plotly")
if (!requireNamespace("mgcv",     quietly = TRUE)) install.packages("mgcv")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
library(plotly)
library(mgcv)
library(reshape2)

# --- 0) Analysis options ----------------------------------------------------
# Target row. A usable target must satisfy FOUR things, and they conflict:
#   1. the package's optimal pair is not (1,2), or views B and C are the same search;
#   2. the B and C PLANES are materially apart -- checked in section 11 by the angle
#      between their normals, NOT by the pair label. Rows 90, 91, 102, 132 and 161
#      all have ~45 deg tilts but their planes coincide EXACTLY (0.00 deg);
#   3. best_pair is reproducible run to run (row 102 flips between (1,2) and (1,3)
#      across identical runs -- see the note's section 4.1);
#   4. the A-B / A-C planes sit near 45 deg, so the section-17 figures are readable.
#
# Criteria (1)-(3) are hard requirements; (4) trades off against how vividly section 18
# can make its point, so two rows are worth keeping in view:
#
#   ROW 6 (set below). Pair (1,3), bit-stable. Its two planes are 90.00 deg apart --
#   orthogonal, the largest separation possible -- so the two slices in section 18 are
#   as different as they can be, separating by 57% of the data's extent. The cost is
#   criterion (4): its A-C plane stands at 88.64 deg, so section 17's A-C figure is
#   near-vertical and awkward to read.
#
#   ROW 104 (kept as the commented alternative). Pair (1,3), planes 18.43 deg apart,
#   tilts 49.18 / 44.74 deg, target 2.694 off the PC1-PC2 plane. Readable section-17
#   figures, but section 18's two slices then separate by only 16% and the contrast is
#   much weaker.
#
# Neither is "correct" -- pick by which section you are looking at. Row 49 fails (1).
# Row 199 fails because the SCRIPT's select_pair() picks (1,2) for it even though the
# package picks (2,3), which collapses views B and C.
i <- 6L #104L

# loan_amnt is the feature dropped, keeping the script-15 / article-section-3 trio.
feat_names <- c("loan_percent_income", "credit_score", "loan_int_rate")

g         <- 60L    # boundary-cloud grid density per axis (g^3 GAM predictions, ONCE)
m_sel     <- 200L   # grid resolution for the pair-selection replication (matches package)
CHECK_PKG <- TRUE   # run the boundarylogic cross-check (sections 8, 9, 12)

# --- 1) Load + filter loan data ---------------------------------------------
loan_raw <- read.csv(
  system.file("extdata", "loan_data.csv", package = "boundarylogic")
)
loan_filtered <- loan_raw[loan_raw$previous_loan_defaults_on_file == "No", ]
loan_filtered <- loan_filtered[, c(feat_names, "loan_status")]
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

# --- 3) Standardize the 3 features ------------------------------------------
sc       <- scale(data[, feat_names], center = TRUE, scale = TRUE)
X_center <- attr(sc, "scaled:center")
X_scale  <- attr(sc, "scaled:scale")
Xstd     <- sc                              # n x 3, standardized
stopifnot(i >= 1L, i <= nrow(Xstd))

# --- 4) PCA (3D -> 3D) via SVD — a pure rotation, nothing discarded ---------
svdx   <- svd(Xstd)
V_full <- svdx$v                            # 3 x 3 orthogonal loading matrix
p      <- ncol(V_full)

# THE PREMISE OF THIS SCRIPT. V is square and orthogonal, so the round trip
# Xstd -> scores -> Xstd is exact for EVERY observation. In 17_2 (4 vars -> 3 PCs)
# this holds only for the rotated target; here it holds for all 221 rows.
Z_full   <- Xstd %*% V_full
rt_resid <- max(abs(Xstd - Z_full %*% t(V_full)))
cat("Round-trip residual (expect ~0):", format(rt_resid, digits = 3), "\n")
stopifnot(
  max(abs(crossprod(V_full) - diag(p))) < 1e-10,   # V orthogonal
  rt_resid < 1e-12                                 # reduction is lossless
)

# --- 5) Fit the nonlinear GAM in standardized 3D space ----------------------
dat_gam <- data.frame(Xstd, class = data$class)
names(dat_gam)[1:3] <- feat_names
gam_fit <- gam(
  class ~ s(loan_percent_income, k = 8) + s(credit_score, k = 8) +
          s(loan_int_rate, k = 8),
  data = dat_gam, family = binomial()
)
summary(gam_fit)

pt_col <- ifelse(dat_gam$class == 1, "red", "blue")   # colour by ACTUAL class

# --- 6) Boundary cloud — computed ONCE, in standardized feature space -------
# The p ~ 0.5 cloud is a set of points in FEATURE space; it does not depend on
# which basis we later draw it in. Computing it once and re-expressing it per
# view (a) avoids running g^3 GAM predictions three times, and (b) guarantees all
# three views show literally the same cloud — nothing is recomputed, so nothing
# can differ between them except the frame.
grid_pca <- expand.grid(
  PC1 = seq(min(Z_full[, 1]), max(Z_full[, 1]), length.out = g),
  PC2 = seq(min(Z_full[, 2]), max(Z_full[, 2]), length.out = g),
  PC3 = seq(min(Z_full[, 3]), max(Z_full[, 3]), length.out = g)
)
Xstd_grid <- as.matrix(grid_pca) %*% t(V_full)   # exact inverse (V is orthogonal)
colnames(Xstd_grid) <- feat_names
grid_p <- predict(gam_fit, newdata = as.data.frame(Xstd_grid), type = "response")

Xstd_cloud <- Xstd_grid[abs(grid_p - 0.5) < 0.02, , drop = FALSE]
cat("Boundary cloud points (p ~ 0.5):", nrow(Xstd_cloud), "\n")
stopifnot(nrow(Xstd_cloud) > 50L)

# --- 7) .bl_rotate_full() — mirror of the package's private .bl_rotate() ----
# The package's .bl_rotate() (R/local_cf.R:45) is the authority for this maths.
# Two deliberate differences here:
#   - it takes the FULL p x p V (rotating a truncated basis would leave the
#     target with a nonzero component outside the drawn subspace);
#   - it returns the FULL rotated basis (p x p), because the 3D plot needs three
#     columns, whereas .bl_rotate() returns only the 2 a biplot needs.
# The algebra is otherwise identical, line for line.
.bl_rotate_full <- function(x_target_st, V, proj_pair) {
  p    <- length(x_target_st)
  Vr   <- V[, proj_pair, drop = FALSE]                  # p x 2

  Y    <- rbind(-x_target_st, rep(0, p), x_target_st)   # 3 x p
  YV   <- Y %*% V                                       # 3 x p
  YVr  <- Y %*% Vr                                      # 3 x 2

  YVr_padded <- cbind(YVr, matrix(0, nrow = 3L, ncol = p - 2L))

  svd_res <- svd(t(YV) %*% YVr_padded)                  # p x p SVD
  A       <- svd_res$v %*% t(svd_res$u)                 # p x p rotation

  V %*% t(A)                                            # rotated full basis (p x p)
}

# Candidate pairs in the package's own order: .generate_pairs() emits
# colexicographically -- for j in 2..p, for i in 1..j-1. For p = 3 that is
# (1,2), (1,3), (2,3). The order matters: it is why max_pairs = 1L pins the
# package to (1,2), which is how view B is cross-checked in section 12.
gen_pairs <- function(p) {
  out <- list()
  for (j in 2L:p) for (ii in 1L:(j - 1L)) out[[length(out) + 1L]] <- c(ii, j)
  do.call(rbind, out)
}
pairs_all  <- gen_pairs(p)
pair_names <- apply(pairs_all, 1L, function(r) sprintf("(%d,%d)", r[1], r[2]))

# --- 8) Package objects (also the source of cutoff / b_margin / metric) -----
# The data handed to bl_prepare_data() is standardised in place below, so the
# package's "raw" feature space IS this script's Xstd space (verified: its
# X_center ~ 0 and X_sd = 1). That makes the replication in section 9 a direct
# comparison rather than a units conversion.
if (CHECK_PKG) {
  bl_dat_std <- bl_dat
  bl_dat_std$train_data[, feat_names] <- Xstd
  bl_dat_std$test_data[,  feat_names] <-
    scale(bl_dat_std$test_data[, feat_names], center = X_center, scale = X_scale)
  bl_dat_std <- bl_set_scaling(bl_dat_std, center = X_center,
                               scale = X_scale, method = "z-score")

  bl_mod <- bl_wrap_model(gam_fit, "GAM", bl_dat_std$var_names,
                          train_data = bl_dat_std$train_data)
  bl_res <- bl_build_result(bl_data = bl_dat_std, bl_model = bl_mod,
                            method = "PCA", title = "Loan status (3 vars)")

  cutoff     <- bl_res$cutoff       # 0.5
  b_margin   <- bl_res$b_margin     # 0.001
  metric_inv <- bl_res$metric_inv
  cat("\nPackage V matches svd(Xstd)$v to:",
      format(max(abs(abs(bl_res$V) - abs(V_full))), digits = 3), "\n")
} else {
  cutoff <- 0.5; b_margin <- 0.001
  # Inverse POOLED WITHIN-CLASS covariance -- the metric the package uses to
  # compare pairs (.compute_metric_inverse: W/(n-g) summed over classes).
  cls_f <- factor(dat_gam$class); W <- matrix(0, p, p)
  for (l in levels(cls_f)) {
    Xk <- Xstd[cls_f == l, , drop = FALSE]
    W  <- W + crossprod(sweep(Xk, 2L, colMeans(Xk)))
  }
  metric_inv <- solve(W / max(nrow(Xstd) - nlevels(cls_f), 1L))
}

# --- 9) How the OPTIMAL PAIR is determined ----------------------------------
# Replicated from bl_find_local_cf(). The rule, in full:
#
#   1. CANDIDATES  .generate_pairs(p, max_pairs = 10L), colex order (section 7).
#   2. ROTATE      .bl_rotate(x_target, V, pair) puts the target IN that pair's
#                  plane and returns Vr_rot (p x 2), tVr_rot (2 x p). Because the
#                  basis is orthonormal, tVr_rot == t(Vr_rot) (checked below).
#   3. GRID        an m x m grid over the padded (10%) range of the rotated
#                  TRAINING scores, mapped back to feature space and predicted.
#   4. CONTOUR     at levels cutoff -/+ b_margin -- TWO levels, not one. Segments
#                  whose implied class equals the target's predicted class are
#                  discarded (a counterfactual must be on the other side), then
#                  actionability filters and the train_ranges test are applied,
#                  then survivors are RE-PREDICTED and kept only if they verify.
#   5. NEAREST POINT WITHIN A PAIR -- plain EUCLIDEAN distance in the rotated 2D
#                  plane (.nearest_idx_block, dx^2 + dy^2). This gives dist_z.
#   6. DISTANCE FOR COMPARING PAIRS -- that point is mapped back to feature space
#                  as B_x, and dist_mah = v' %*% metric_inv %*% v, v = B_x - x_obs,
#                  where metric_inv is the inverse POOLED WITHIN-CLASS covariance.
#                  NOTE: this is the SQUARED Mahalanobis distance -- no sqrt. The
#                  package's own summary prints sqrt() of it.
#   7. SELECTION   argmin over pairs of dist_mah (distance = "mahalanobis", the
#                  default), else of dist_z.
#
# So TWO DIFFERENT METRICS act at TWO DIFFERENT STAGES: Euclidean-in-plane picks
# the point, Mahalanobis-in-feature-space picks the pair. They can disagree --
# for row 49 of this data Euclidean would choose (1,3) while Mahalanobis chooses
# (1,2). This is the most misreadable part of the algorithm.
#
# Replicated here: steps 1-3, 5-7, and the contour/segment-class part of 4. NOT
# replicated: the actionability filters (none are set) and the train_ranges test.
#
# WARNING -- the train_ranges omission is NOT harmless. An earlier version of this
# comment called it "a no-op for an in-sample target". That was wrong: the test
# filters the BOUNDARY CANDIDATES, not the target, so whether the target is
# in-sample is irrelevant. train_ranges is the observed min/max box of the
# standardised training data, and the search grid deliberately extends 10% beyond
# the score range, so mapped-back candidates routinely fall outside it. Measured
# fraction of candidates removed:
#
#             pair (1,2)   (1,3)   (2,3)
#   row   6       31%      18%     56%
#   row 199       56%      51%     49%
#
# Dropping up to 56% of the candidates changes which one is nearest, hence the
# distances, hence the selected pair. That is why the comparison table in section
# 15 usually DISAGREES with the package -- roughly 3 rows in 13 agree. Read that
# table as a MEASUREMENT OF THIS OMISSION, not as a fidelity check. (Rotating from
# bl_res$V rather than svd(Xstd)$v does not help: 2/13. The V source is not the
# cause.) Nothing drawn depends on select_pair() any more -- views B and C come
# from the package's own Vr_rot in section 11 -- so this stays an exhibit.
x_t <- as.numeric(Xstd[i, ])

select_pair <- function() {
  d_z   <- setNames(rep(NA_real_, nrow(pairs_all)), pair_names)
  d_mah <- d_z

  x_t_df <- as.data.frame(matrix(x_t, nrow = 1L)); colnames(x_t_df) <- feat_names
  tgt_p     <- as.numeric(predict(gam_fit, newdata = x_t_df, type = "response"))
  tgt_class <- as.integer(tgt_p >= cutoff)
  cat(sprintf("Target row %d: p = %.4f -> class %d\n", i, tgt_p, tgt_class))

  for (k in seq_len(nrow(pairs_all))) {
    pr     <- pairs_all[k, ]
    Vrho   <- .bl_rotate_full(x_t, V_full, pr)
    Vr_rot <- Vrho[, 1:2, drop = FALSE]

    Z_target <- as.numeric(x_t %*% Vr_rot)
    Z_train  <- Xstd %*% Vr_rot
    zr   <- range(Z_train); zpad <- diff(zr) * 0.1
    sq   <- seq(zr[1] - zpad, zr[2] + zpad, length.out = m_sel)
    Zg   <- as.matrix(expand.grid(sq, sq))
    Xg   <- as.data.frame(Zg %*% t(Vr_rot)); colnames(Xg) <- feat_names
    pg   <- predict(gam_fit, newdata = Xg, type = "response")
    ct   <- contourLines(sq, sq, matrix(pg, ncol = m_sel, byrow = FALSE),
                         levels = c(cutoff - b_margin, cutoff + b_margin))
    if (length(ct) == 0L) next

    keep <- list()
    for (seg in ct) {
      if (as.integer(seg$level >= cutoff) == tgt_class) next   # wrong side
      keep[[length(keep) + 1L]] <- cbind(seg$x, seg$y)
    }
    if (length(keep) == 0L) next
    Mi <- do.call(rbind, keep)

    d2  <- (Mi[, 1] - Z_target[1])^2 + (Mi[, 2] - Z_target[2])^2   # EUCLIDEAN, in plane
    idx <- which.min(d2)
    d_z[k] <- sqrt(d2[idx])

    B_x      <- as.numeric(Mi[idx, , drop = FALSE] %*% t(Vr_rot))  # back to feature space
    v        <- B_x - x_t
    d_mah[k] <- as.numeric(t(v) %*% metric_inv %*% v)              # SQUARED Mahalanobis
  }
  list(d_z = d_z, d_mah = d_mah, best = pairs_all[which.min(d_mah), ])
}

sel       <- select_pair()
best_pair <- as.integer(sel$best)
cat(sprintf("Script's optimal pair: (%d,%d)\n", best_pair[1], best_pair[2]))

if (identical(best_pair, c(1L, 2L)))
  warning(sprintf(
    paste("Target row %d's optimal pair IS (1,2), so views B and C are the SAME",
          "rotation and the comparison shows nothing. Pick a different row --",
          "row 6 gives (1,3) on this data."), i), call. = FALSE)

# --- 10) make_view() — one code path, three bases ---------------------------
# The comparison is only trustworthy if the views differ in NOTHING but the
# basis, so all three are built here rather than in three copy-pasted blocks.
make_view <- function(V_draw, view_title, axis_lab, proj_pair = NULL) {
  Z  <- as.data.frame(Xstd %*% V_draw);        names(Z)  <- c("s1", "s2", "s3")
  Cl <- as.data.frame(Xstd_cloud %*% V_draw);  names(Cl) <- c("s1", "s2", "s3")

  surf_fit <- gam(s3 ~ s(s1, s2, k = 60), data = Cl)
  s1_seq <- seq(min(Z$s1), max(Z$s1), length.out = 80)
  s2_seq <- seq(min(Z$s2), max(Z$s2), length.out = 80)
  gr <- expand.grid(s1 = s1_seq, s2 = s2_seq)
  gr$s3 <- predict(surf_fit, newdata = gr)

  # UNMASKED matrix for the contour (masking first would truncate the curve at
  # the range edges); masked copy for display. acast() returns s2 x s1 but
  # contourLines() wants z indexed [x, y] = [s1, s2] -- hence t().
  surf_full <- acast(gr, s2 ~ s1, value.var = "s3")
  gr$s3[gr$s3 < min(Z$s3) | gr$s3 > max(Z$s3)] <- NA
  surf_disp <- acast(gr, s2 ~ s1, value.var = "s3")

  cls <- contourLines(x = s1_seq, y = s2_seq, z = t(surf_full), levels = 0)

  pc_t    <- as.numeric(Z[i, ])
  ss_keep <- sum((Xstd %*% V_draw)^2) / sum(Xstd^2)
  ang     <- if (is.null(proj_pair)) NA_real_ else {
    u <- as.numeric(x_t %*% V_full)
    acos(min(1, sqrt(sum(u[proj_pair]^2)) / sqrt(sum(u^2)))) * 180 / pi
  }

  cat(sprintf(
    "%-48s | off-plane |s3| = %9.3e | SS retained = %.4f | angle = %s | segments = %d\n",
    view_title, abs(pc_t[3]), ss_keep,
    if (is.na(ang)) "    n/a" else sprintf("%5.1f deg", ang), length(cls)))

  arrows <- function(fig, twoD) {
    for (k in seq_len(nrow(V_draw))) {
      if (twoD) {
        fig <- fig |> add_trace(
          x = c(0, V_draw[k, 1]) * 3, y = c(0, V_draw[k, 2]) * 3,
          type = "scatter", mode = "lines+text",
          text = c("", feat_names[k]), textposition = "top center",
          textfont = list(size = 14, color = "#222222"),
          line = list(color = "#5E88FC", width = 5), showlegend = FALSE)
      } else {
        fig <- fig |> add_trace(
          x = c(0, V_draw[k, 1]) * 3, y = c(0, V_draw[k, 2]) * 3,
          z = c(0, V_draw[k, 3]) * 3,
          type = "scatter3d", mode = "lines+text",
          text = c("", feat_names[k]), textposition = "top center",
          textfont = list(size = 14, color = "#222222"),
          line = list(color = "#5E88FC", width = 5), showlegend = FALSE)
      }
    }
    fig
  }

  # ---- 3D ----
  f3 <- plot_ly() |>
    add_markers(data = Z, x = ~s1, y = ~s2, z = ~s3,
                marker = list(size = 5, color = pt_col, opacity = 0.4),
                name = "Data points", showlegend = FALSE) |>
    add_surface(x = s1_seq, y = s2_seq, z = surf_disp,
                colorscale = list(c(0, "green"), c(1, "green")),
                opacity = 0.4, showscale = FALSE,
                name = "GAM decision surface (p = 0.5)")
  f3 <- arrows(f3, twoD = FALSE)
  for (cl in cls) {
    f3 <- f3 |> add_trace(x = cl$x, y = cl$y, z = rep(0, length(cl$x)),
                          type = "scatter3d", mode = "lines",
                          line = list(width = 6, color = "black"),   # in-plane decision boundary
                          name = "Boundary at third axis = 0", showlegend = FALSE)
  }
  f3 <- f3 |>
    add_markers(x = pc_t[1], y = pc_t[2], z = pc_t[3],
                marker = list(size = 9, symbol = "diamond", color = "magenta"),
                name = "Target", showlegend = FALSE) |>
    layout(title = paste0(view_title, " — 3D"),
           scene = list(aspectmode = "data",
                        xaxis = list(title = axis_lab[1]),
                        yaxis = list(title = axis_lab[2]),
                        zaxis = list(title = axis_lab[3]),
                        camera = list(eye = list(x = 1.4, y = 1.4, z = 0.9))),
           showlegend = FALSE)

  # ---- 2D biplot: the same object viewed down the third axis ----
  f2 <- plot_ly() |>
    add_markers(data = Z, x = ~s1, y = ~s2,
                marker = list(size = 7, color = pt_col, opacity = 0.4),
                name = "Data points", showlegend = FALSE)
  f2 <- arrows(f2, twoD = TRUE)
  for (cl in cls) {
    f2 <- f2 |> add_trace(x = cl$x, y = cl$y, type = "scatter", mode = "lines",
                          line = list(width = 6, color = "black"),   # in-plane decision boundary
                          name = "Decision boundary", showlegend = FALSE)
  }
  f2 <- f2 |>
    add_markers(x = pc_t[1], y = pc_t[2],
                marker = list(size = 13, symbol = "diamond", color = "magenta"),
                name = "Target", showlegend = FALSE) |>
    layout(title = paste0(view_title, " — 2D biplot"),
           xaxis = list(title = axis_lab[1], scaleanchor = "y", scaleratio = 1),
           yaxis = list(title = axis_lab[2]),
           showlegend = FALSE)

  list(fig3d = f3, fig2d = f2, Z = Z, cls = cls, pc_t = pc_t,
       V = V_draw, ss = ss_keep)
}

# --- 11) The three views ----------------------------------------------------
cat("\n--- Views (SS retained is 1.0000 for all three: the reduction is lossless) ---\n")

# WHERE THE DRAWN PLANES COME FROM. Views B and C are built from the PACKAGE's own
# rotated bases (bl_loc$Vr_rot), not from this script's select_pair(). Two reasons:
#   - select_pair() omits the train_ranges filter (section 9), so its chosen pair
#     disagrees with the package on most rows -- the figures should not inherit that;
#   - Vr_rot IS the basis the package's local biplot uses, so the views are then the
#     real thing rather than a reconstruction of it. Script 17_2 already works this way.
# Vr_rot is p x 2; the 3D view needs a third column, so take the direction in the
# orthogonal complement that retains the most variance (lifted from 17_2 section 7).
third_axis <- function(Vr) {
  R  <- Xstd - Xstd %*% Vr %*% t(Vr)
  sv <- svd(R)
  stopifnot(sv$d[1] > 1e-8)      # degenerate only if the data lie in the plane
  sv$v[, 1]
}

V_A <- V_full
if (CHECK_PKG) {
  tgt    <- bl_select_target(bl_res, target = bl_res$train_data[i, ])
  loc_12 <- bl_find_local_cf(bl_res, tgt, max_pairs = 1L, verbose = FALSE)  # pinned to (1,2)
  loc_op <- bl_find_local_cf(bl_res, tgt, verbose = FALSE)                  # full search
  stopifnot(identical(as.integer(loc_12$best_pair), c(1L, 2L)),
            isTRUE(loc_12$solution_found), isTRUE(loc_op$solution_found))
  pkg_pair <- as.integer(loc_op$best_pair)
  V_B <- cbind(loc_12$Vr_rot, third_axis(loc_12$Vr_rot))
  V_C <- cbind(loc_op$Vr_rot, third_axis(loc_op$Vr_rot))
  cat(sprintf("Planes drawn from the package: B = pinned (1,2), C = optimal (%d,%d)\n",
              pkg_pair[1], pkg_pair[2]))
} else {
  # No package objects available -- fall back to the script's own rotation. The
  # views are then this script's planes, NOT the package's; section 15c shows the
  # two can differ substantially.
  pkg_pair <- best_pair
  V_B <- .bl_rotate_full(x_t, V_full, c(1L, 2L))
  V_C <- .bl_rotate_full(x_t, V_full, best_pair)
  cat("CHECK_PKG is FALSE: planes are the SCRIPT's own, not the package's.\n")
}
stopifnot(max(abs(crossprod(V_B) - diag(3))) < 1e-8,
          max(abs(crossprod(V_C) - diag(3))) < 1e-8)

# Axis labels: view A is genuinely unrotated, so PC1/PC2/PC3 are correct there.
# The rotated views get RC1/RC2/RC3 -- see the header. The named pair is real
# provenance for WHICH PLANE was rotated, so it goes in the title, not on an axis.
view_A <- make_view(V_A, "A. Unrotated PCA biplot", c("PC1", "PC2", "PC3"))
view_B <- make_view(V_B,
                    "B. RC biplot - rotated from the (1,2) plane",
                    c("RC1", "RC2", "RC3 (residual)"), c(1L, 2L))
view_C <- make_view(V_C,
                    sprintf("C. RC biplot - rotated from the (%d,%d) plane [optimal]",
                            pkg_pair[1], pkg_pair[2]),
                    c("RC1", "RC2", "RC3 (residual)"), pkg_pair)

# COLLAPSE GUARD -- on the PLANE, not the pair label. The label test in section 9
# cannot see the real failure: rows 137, 90, 91, 102, 132 and 161 all carry a pair
# label other than (1,2) while their B and C planes coincide exactly, so the two
# views are the same picture drawn twice.
nB_chk <- as.numeric(t(V_A) %*% V_B[, 3])
nC_chk <- as.numeric(t(V_A) %*% V_C[, 3])
bc_sep <- acos(min(1, abs(sum(nB_chk * nC_chk)))) * 180 / pi
cat(sprintf("B and C plane separation: %.2f deg\n", bc_sep))
if (bc_sep < 5)
  warning(sprintf(paste("B and C planes are only %.2f deg apart -- the two views are",
                        "effectively the same picture. Pick a different target row."),
                  bc_sep), call. = FALSE)

# The rotation's whole job: put the target ON the drawn plane. It is off-plane in
# A (that is the motivation) and exactly on-plane in B and C.
stopifnot(
  abs(view_A$pc_t[3]) > 1e-6,   # A: target floats off the plane
  abs(view_B$pc_t[3]) < 1e-8,   # B: on the plane
  abs(view_C$pc_t[3]) < 1e-8,   # C: on the plane
  abs(view_A$ss - 1) < 1e-10, abs(view_B$ss - 1) < 1e-10, abs(view_C$ss - 1) < 1e-10,
  length(view_A$cls) >= 1L, length(view_B$cls) >= 1L, length(view_C$cls) >= 1L
)

# --- 12) The rotation itself, and what the pair label really means ----------
# The rotation is pinned on exactly ONE degree of freedom (see section 15c). The
# constraint is A'u ~ w, where u = V'x is the target in PC coordinates and w is u
# with the non-pair entries zeroed. Since <u,w> = ||w||^2, cos(theta) = ||w||/||u||.
# A is orthogonal, so A'u = ||u|| * w/||w||: the target swings along an ARC OF
# CONSTANT RADIUS ||u||.
#
# So the rotation acts entirely inside span(u, w) -- a single 2D plane -- and that
# is what the panel below draws, in the orthonormal frame e1 = w/||w||, e2 = the
# part of u orthogonal to w. This is A3's panel C generalised from p = 2 to any p.
#
# The point it makes visually: ORTHOGONAL PROJECTION (what the unrotated biplot
# does) drops the target onto the axis at ||w||, SHORTENING it. The ROTATION instead
# brings the plane up to the target, landing it at ||u|| -- full length preserved,
# which is the x_hat = x identity.
rotation_panel <- function(pair, title) {
  u  <- as.numeric(x_t %*% V_full)          # target in PC coordinates
  w  <- u; w[-pair] <- 0                    # its shadow on the named plane
  un <- sqrt(sum(u^2)); wn <- sqrt(sum(w^2))
  th <- acos(min(1, wn / un))               # the rotation angle
  h  <- sqrt(max(0, un^2 - wn^2))           # height of the target above the plane

  a   <- if (th > 1e-9) seq(0, th, length.out = 200L) else c(0, 0)
  lim <- un * 1.18

  plot_ly() |>
    # the named plane, seen edge-on in this frame
    add_trace(x = c(-lim, lim), y = c(0, 0), type = "scatter", mode = "lines",
              line = list(color = "#bbbbbb", width = 2), showlegend = FALSE) |>
    # the constant-radius arc the target travels along
    add_trace(x = un * cos(a), y = un * sin(a), type = "scatter", mode = "lines",
              line = list(color = "#888888", width = 3, dash = "dot"),
              showlegend = FALSE) |>
    # the target's vector before rotation
    add_trace(x = c(0, wn), y = c(0, h), type = "scatter", mode = "lines",
              line = list(color = "magenta", width = 3), showlegend = FALSE) |>
    # the orthogonal-projection drop -- what the UNROTATED biplot does
    add_trace(x = c(wn, wn), y = c(h, 0), type = "scatter", mode = "lines",
              line = list(color = "#888888", width = 2, dash = "dash"),
              showlegend = FALSE) |>
    add_trace(x = wn, y = 0, type = "scatter", mode = "markers+text",
              marker = list(size = 11, symbol = "square", color = "#888888"),
              text = sprintf("orthogonal projection (%.3f)", wn),
              textposition = "bottom center",
              textfont = list(size = 13, color = "#555555"), showlegend = FALSE) |>
    # target before (off the plane) and after (on it, at full length)
    add_trace(x = wn, y = h, type = "scatter", mode = "markers+text",
              marker = list(size = 13, symbol = "diamond", color = "magenta"),
              text = "target (before)", textposition = "top center",
              textfont = list(size = 13, color = "magenta"), showlegend = FALSE) |>
    add_trace(x = un, y = 0, type = "scatter", mode = "markers+text",
              marker = list(size = 13, symbol = "circle", color = "magenta"),
              text = sprintf("target after rotation (%.3f)", un),
              textposition = "top center",
              textfont = list(size = 13, color = "magenta"), showlegend = FALSE) |>
    layout(
      title = sprintf("%s  -  theta = %.1f deg", title, th * 180 / pi),
      xaxis = list(title = "e1 = w / ||w||  (inside the named plane)",
                   scaleanchor = "y", scaleratio = 1, zeroline = FALSE),
      yaxis = list(title = "e2  (orthogonal to the plane)", zeroline = FALSE),
      showlegend = FALSE)
}

# Evidence that the named pair is NOT an axis identity. For each candidate pair:
# theta, the principal angles between the rotated plane and the eigenplane the pair
# NAMES, and how much of that plane's second direction actually lies inside it.
prin_ang <- function(Q1, Q2) {
  sv <- svd(t(Q1) %*% Q2)$d
  acos(pmin(1, pmax(-1, sv))) * 180 / pi
}
cat("\n--- 12) What the pair label is worth ---\n")
pl <- data.frame(pair = pair_names, theta_deg = NA_real_,
                 ang1_deg = NA_real_, ang2_deg = NA_real_,
                 pct_in_named = NA_real_, row.names = NULL)
for (k in seq_len(nrow(pairs_all))) {
  pr <- pairs_all[k, ]
  Vd <- .bl_rotate_full(x_t, V_full, pr)
  uu <- as.numeric(x_t %*% V_full); ww <- uu; ww[-pr] <- 0
  pl$theta_deg[k] <- acos(min(1, sqrt(sum(ww^2)) / sqrt(sum(uu^2)))) * 180 / pi
  ang <- prin_ang(Vd[, 1:2, drop = FALSE], V_full[, pr, drop = FALSE])
  pl$ang1_deg[k] <- ang[1]; pl$ang2_deg[k] <- ang[2]
  # the plane's 2nd direction, taken orthogonal to the target direction
  tv <- as.numeric(x_t %*% Vd[, 1:2]); td <- Vd[, 1:2] %*% (tv / sqrt(sum(tv^2)))
  d2 <- Vd[, 1:2] %*% c(0, 1); d2 <- d2 - td * sum(d2 * td)
  d2 <- d2 / sqrt(sum(d2^2))
  Pn <- V_full[, pr, drop = FALSE] %*% t(V_full[, pr, drop = FALSE])
  pl$pct_in_named[k] <- 100 * sum((Pn %*% d2)^2)
  # which original PC is each rotated column actually closest to?
  cl <- function(col) which.min(acos(pmin(1, abs(as.numeric(t(col) %*% V_full)))))
  cat(sprintf("  pair %s: RC1 is closest to PC%d, RC2 to PC%d\n",
              pair_names[k], cl(Vd[, 1]), cl(Vd[, 2])))
}
print(pl, digits = 4)
# CAREFUL WITH ang1 HERE. In 3D it is 0 for EVERY pair because any two 2-planes in
# 3-space must intersect in a line -- it is forced by the dimension, not produced by
# the rotation (two RANDOM planes give 0.00, 83.34 / 0.00, 47.79 / 0.00, 87.81). Nor
# is the shared line the target's direction: only 0.4435 of that direction lies in
# the PC1-PC2 plane, which is exactly its unrotated sample predictivity in section 13.
# CORRECTED 2026-09-10. An earlier version of this comment said the note's "the
# first principal angle is always 0 -- the target's own direction is always shared"
# was "a statement about p = 4, where a 0 really does mean the target direction".
# That conflated two DIFFERENT comparisons, and the p = 4 measurement settles it:
#   THIS TABLE compares a rotated plane against the eigenplane its label NAMES.
#     That eigenplane does not contain the target, so nothing forces a zero -- and
#     at p = 4 there is none: ang1 runs 1.11 to 55.62 deg across the six pairs
#     (17_2 section 9, 17_3 section 2). Here at p = 3 it is 0 only because any two
#     2-planes in 3-space must share a line.
#   THE NOTE'S SECTION 4.1 TABLE compares perturbed against unperturbed rotations
#     of ONE pair. Both of those contain the target, so ang1 = 0 necessarily, and
#     at p = 4 that zero IS meaningful. See 17_2 section 8b.
# So the zero is meaningful for that comparison, not for this one, at either p.
cat("  ang1 is 0 for every pair, but in 3D that is forced (any two planes share a\n")
cat("  line) and the shared line is NOT the target direction -- see the comment above.\n")
cat("  ang2 is the informative one.\n")
# theta does NOT predict pct_in_named. Measured across all 663 (row, pair)
# combinations on this data: Pearson r = +0.33, Spearman +0.14, and among the 74
# cases with theta < 5 deg the overlap spans the FULL range 0% to 100%. Row 6's
# tidy pairing of theta = 3.03 deg with 99.9% overlap was a coincidence; row 137
# has theta = 1.49 deg with 0.1% overlap. The completion is simply arbitrary, so
# the overlap has to be read off the table -- it cannot be inferred from theta.
cat("  theta does NOT predict pct_in_named (r = +0.33; among theta < 5 deg the\n")
cat("  overlap spans 0-100%). Read the overlap off the table, do not infer it.\n")

# --- 13) What the rotation costs -------------------------------------------
# A3's diagnostic block (its section 4), with ONE correction: RSS must use the
# TRACE. A3 as written summed EVERY ENTRY of X'X and of X'X V J V', which is not a
# sum of squares and made the rotation appear to IMPROVE global fit -- the opposite
# of the method's actual trade-off. Corrected in A3 too (2026-09-03). With the
# trace, tr(X'X) = (n-1)*p for standardised data, and the unrotated RSS equals the
# dropped component's variance exactly; both are asserted below.
#
# J = diag(1,1,0) selects the 2D biplot: keep the first two axes, drop the third.
D <- svdx$d; S <- D^2
stopifnot(min(D) > 1e-8)
J      <- diag(c(1, 1, 0))
XtX    <- crossprod(Xstd)
tot_SS <- sum(diag(XtX))

rc_diagnostics <- function(V_draw, label) {
  # Uv uses the ORIGINAL D. That is fine: J is diagonal so it commutes with D and
  # the D factors cancel, leaving ||proj of x_i onto the kept plane||^2 / ||x_i||^2
  # -- genuine sample predictivity. The target hitting exactly 1 confirms it.
  Uv  <- Xstd %*% V_draw %*% solve(diag(D))
  adq <- diag(V_draw %*% J %*% t(V_draw)) / diag(V_draw %*% t(V_draw))
  axp <- diag(V_draw %*% diag(S) %*% J %*% t(V_draw)) /
         diag(V_draw %*% diag(S) %*% t(V_draw))
  smp <- diag(Uv %*% diag(S) %*% J %*% t(Uv)) / diag(Uv %*% diag(S) %*% t(Uv))
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

cat("\n--- 13) The cost of the rotation (J = diag(1,1,0): the 2D biplot) ---")
dg_A <- rc_diagnostics(V_A, "A. Unrotated PCA")
dg_B <- rc_diagnostics(V_B, "B. RC biplot, rotated from (1,2)")
dg_C <- rc_diagnostics(V_C, sprintf("C. RC biplot, rotated from (%d,%d) [optimal]",
                                    best_pair[1], best_pair[2]))

# NOTE: `adequacy` indexes the ROWS of V, so it is PER VARIABLE (per feature), not
# per axis. An earlier version of this line read "the third axis stops meaning PC3",
# which mislabelled one variable's representation as a property of an axis. The
# effect is a REDISTRIBUTION across variables -- some improve, some collapse -- so
# the collapsing variable is computed rather than assumed to be the third.
#
# Compare the two rotated rows printed above: they do not sacrifice the same variable.
# Which one collapses depends on where the rotated plane ends up, not on any ranking
# of the variables.
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
  abs(dg_A$rss - S[3]) < 1e-8,      # unrotated RSS IS the dropped component's variance
  # The unrotated PC1-PC2 plane minimises the residual over all rank-2 planes, so any
  # rotation of it must increase SS. ("optimal" is reserved in this script for
  # bl_find_local_cf()'s selection, so it is not used for that least-squares sense.)
  dg_B$rss > dg_A$rss,
  dg_C$rss > dg_A$rss,
  abs(dg_B$smp[i] - 1) < 1e-8,      # the defining property: the target becomes exact
  abs(dg_C$smp[i] - 1) < 1e-8,
  dg_B$xh < 1e-8, dg_C$xh < 1e-8,   # x_hat = x for the target ...
  dg_A$xh > 1e-3                    # ... and NOT before the rotation
)

# --- 14) The 3D and the 2D view are the same object ------------------------
# The 2D biplot IS the 3D figure looked at down the third axis. What makes that
# reading EXACT for the target -- and only for the target -- is that the rotation
# puts it at RC3 = 0. Every other point keeps a nonzero RC3, so for them the 2D
# panel is a projection, and their sample predictivity in section 13 says by how
# much. This is the r = 2 form of "PC2 disappears" in the 1D A3 example.
cat("\n--- 14) 3D vs 2D: exact for the target, a projection for everyone else ---\n")
for (nm in c("A", "B", "C")) {
  vw <- switch(nm, A = view_A, B = view_B, C = view_C)
  cat(sprintf("  view %s: target third-axis coord = %9.3e | points within 1e-8 of the plane: %d of %d\n",
              nm, vw$pc_t[3], sum(abs(vw$Z$s3) < 1e-8), nrow(vw$Z)))
}
stopifnot(
  # exactly one point lies in the drawn plane after rotation: the target
  sum(abs(view_B$Z$s3) < 1e-8) == 1L,
  sum(abs(view_C$Z$s3) < 1e-8) == 1L,
  # and it IS the target, not some other point that happens to sit there
  which.min(abs(view_C$Z$s3)) == i
)
cat(sprintf(paste0("  -> looking down RC3 reproduces the 2D biplot EXACTLY for the target,\n",
                   "     and as an orthogonal projection for the other %d points.\n"),
            nrow(view_C$Z) - 1L))

# --- 15) Cross-check against boundarylogic ----------------------------------
# Read section 15c before trusting any plane comparison: the rotated plane is
# NOT uniquely determined, so "does my plane equal the package's plane?" is the
# wrong question to assert on. 15a asks the right one.
if (CHECK_PKG) {
  # tgt / loc_12 / loc_op were created in section 11, where they define the drawn
  # planes. Reused here rather than recomputed, so the comparison is against the
  # exact objects the figures use.
  P <- function(M) M %*% t(M)   # plane projector: invariant to in-plane basis + sign

  # ---- 15-0) Pair selection: the script's replication vs the package -------
  cat("
--- Pair selection: script vs package ---
")
  cmp <- data.frame(
    pair          = pair_names,
    script_d_z    = round(sel$d_z, 6),
    pkg_d_z       = round(as.numeric(loc_op$all_distances[pair_names]), 6),
    script_d_mah2 = round(sel$d_mah, 6),
    pkg_d_mah2    = round(as.numeric(loc_op$all_distances_mahalanobis[pair_names]), 6),
    row.names     = NULL
  )
  cmp$selected <- ifelse(cmp$pair == sprintf("(%d,%d)", best_pair[1], best_pair[2]),
                         "<-", "")
  print(cmp)
  cat("(d_mah2 is the SQUARED Mahalanobis distance -- the package stores it unrooted.)
")
  # Rows agree to ~0.3% where the plane happens to coincide, and diverge where it
  # does not -- pair (1,2) here. That divergence is the section 15c effect, not a
  # replication error. The residual gap on the matching rows is the two filters
  # this script deliberately omits (actionability, train_ranges).
  if (!identical(as.integer(loc_op$best_pair), best_pair)) {
    warning(sprintf(
      "Script chose (%d,%d) but the package chose (%d,%d).",
      best_pair[1], best_pair[2],
      loc_op$best_pair[1], loc_op$best_pair[2]), call. = FALSE)
  } else {
    cat("Pair selection AGREES with the package.
")
  }

  # ---- 15a) Is .bl_rotate_full() a faithful mirror of .bl_rotate()? --------
  # THE meaningful test: same x_target, same V, same pair -> same plane. Any
  # comparison against the package's LIVE objects instead varies its inputs too,
  # which 15c shows is not a fair test.
  cat("
--- 15a) Mirror test vs .bl_rotate(), identical inputs ---
")
  mirror_max <- 0
  for (k in seq_len(nrow(pairs_all))) {
    pr   <- pairs_all[k, ]
    mine <- .bl_rotate_full(x_t, V_full, pr)[, 1:2, drop = FALSE]
    pkgr <- boundarylogic:::.bl_rotate(x_t, V_full, pr)$Vr_rot
    dd   <- max(abs(P(mine) - P(pkgr)))
    mirror_max <- max(mirror_max, dd)
    cat(sprintf("  pair %s plane diff: %.3e
", pair_names[k], dd))
  }
  stopifnot(mirror_max < 1e-12)
  cat("  .bl_rotate_full() reproduces .bl_rotate() exactly.
")

  # ---- 15b) What the rotation actually guarantees --------------------------
  # Not "a particular plane", but: the target lies IN the plane, so its two
  # in-plane coordinates carry its full length. This is the invariant that makes
  # the local biplot exact for the target -- and it holds for every pair.
  cat("
--- 15b) Invariants (||x_t|| =", format(sqrt(sum(x_t^2)), digits = 8), ") ---
")
  for (nm in c("B", "C")) {
    Vd <- if (nm == "B") V_B else V_C
    tp <- as.numeric(x_t %*% Vd[, 1:2, drop = FALSE])
    cat(sprintf("  view %s: in-plane norm = %.8f | off-plane |s3| = %.3e | orthonormal: %s
",
                nm, sqrt(sum(tp^2)), abs(as.numeric(x_t %*% Vd)[3]),
                max(abs(crossprod(Vd) - diag(p))) < 1e-10))
    stopifnot(abs(sqrt(sum(tp^2)) - sqrt(sum(x_t^2))) < 1e-8)
  }
  cat("  Target length is fully carried by the plane in both -> it lies exactly in it.
")

  # ---- 15c) FINDING: the rotated plane is NOT unique -----------------------
  # .bl_rotate solves an orthogonal Procrustes problem whose matrix is
  # t(YV) %*% YVr_padded. Because Y = rbind(-x, 0, x) has rank 1, that matrix is
  # RANK 1 (singular values below): A = svd$v %*% t(svd$u) is therefore pinned on
  # only ONE direction, and the remaining p-1 are an arbitrary orthonormal
  # completion. Consequence: any plane containing the target is a valid answer.
  # The completion is deterministic for identical inputs (12a) but NOT stable to
  # perturbation -- and bl_res$V differs from svd(Xstd)$v at ~1e-15, which is
  # enough to land on a different plane. It did for (1,2) and did not for (1,3).
  #
  # This matters beyond aesthetics: the nearest boundary point depends on the
  # plane, so dist_z and the squared Mahalanobis distance -- and hence pair
  # selection -- inherit this sensitivity. See the (1,2) row of the table above.
  Yq   <- rbind(-x_t, rep(0, p), x_t)
  Mq   <- t(Yq %*% V_full) %*% cbind(Yq %*% V_full[, 1:2], matrix(0, 3L, p - 2L))
  cat("
--- 15c) Procrustes matrix singular values (pair (1,2)):",
      format(svd(Mq)$d, digits = 3), "---
")
  cat("  numerical rank:", sum(svd(Mq)$d > 1e-10 * max(svd(Mq)$d)), "of", p,
      "-> the plane is free to rotate about the target direction.
")
  cat(sprintf("  plane vs package  B (1,2): %.3e   C (%d,%d): %.3e
",
              max(abs(P(V_B[, 1:2, drop = FALSE]) - P(loc_12$Vr_rot))),
              best_pair[1], best_pair[2],
              max(abs(P(V_C[, 1:2, drop = FALSE]) - P(loc_op$Vr_rot)))))
  cat("  (Not asserted: a nonzero value here is the free completion, not an error.)
")
  cat("  tVr_rot == t(Vr_rot)?",
      max(abs(loc_op$tVr_rot - t(loc_op$Vr_rot))) < 1e-12, "
")
}

# --- 16) Figures ------------------------------------------------------------
# THE ROTATION ITSELF (section 12) -- the A3 panel-C analogue. Run these first:
# the (1,2) panel shows a large, obvious rotation, the optimal-pair panel a nearly
# degenerate one, and that contrast is why the pair label is worth more in one case
# than the other.
rot_panel_B <- rotation_panel(c(1L, 2L), "Rotation onto the (1,2) plane")
rot_panel_C <- rotation_panel(best_pair,
                 sprintf("Rotation onto the (%d,%d) plane [optimal]",
                         best_pair[1], best_pair[2]))
rot_panel_B; rot_panel_C

# Script's own views (run these one at a time in RStudio). For each rotated view,
# the 3D figure viewed down RC3 IS the 2D biplot beside it (section 14).
view_A$fig3d; view_A$fig2d
view_B$fig3d; view_B$fig2d
view_C$fig3d; view_C$fig2d

# The package's own biplots, for side-by-side comparison with view C's 2D panel:
if (CHECK_PKG) {
  plot(bl_res)     # global PCA biplot with the decision boundary
  plot(loc_12)
  plot(loc_op)     # rotated LOCAL biplot: boundary + nearest counterfactual
}

# --- 17) "A-B plane" and "A-C plane": B and C shown INSIDE the unrotated view
# ----------------------------------------------------------------------------
# Views A, B and C are three separate figures, so the relationship between them has
# to be taken on trust: nothing in A shows WHERE B's or C's biplot plane actually
# sits. These two figures fix that. Each is view A's 3D plot with the other view's
# biplot plane drawn into it -- tilt the plot until that plane is horizontal and you
# are looking at B (or C).
#
# THE CONSTRUCTION IS EXACT, NOT A LIKENESS.
# View A plots a = x'V_A; view B plots b = x'V_B. V_A is orthogonal, so x = V_A a and
#
#     b = a' (V_A' V_B)   =>   e1 = V_A'V_B[,1],  e2 = V_A'V_B[,2],  n = V_A'V_B[,3]
#
# n is B's residual axis expressed in A's coordinates, so B's biplot plane seen from
# inside view A is the plane through the ORIGIN with normal n and in-plane directions
# e1, e2. [e1 e2 n] is orthonormal, and A's coordinates multiplied by it reproduce
# view B's coordinates EXACTLY (asserted below, ~1e-15). Rotating A by that matrix
# does not merely resemble B -- it IS B.
#
# TWO THINGS THE FIGURE SHOWS THAT THE SEPARATE VIEWS CANNOT.
#   1. The plane passes exactly through the TARGET. That is the rotation's defining
#      property (x_t'V_B[,3] = 0) seen from the other side.
#   2. How far the plot must be TIPPED to reach B or C -- ~86 deg and ~89 deg here.
#      That angle is the visual measure of how badly the unrotated biplot
#      misrepresents this target.
#
# WHY THE IN-PLANE AXES ARE DRAWN. Making the plane horizontal fixes only two of the
# three rotational degrees of freedom; a spin about the vertical axis stays free. The
# e1/e2 lines pin it, and they are exactly B's (or C's) RC1 and RC2 directions, so the
# correspondence is exact rather than up-to-a-spin.
#
# NOTE ON THE TILT ANGLES vs SECTION 12. The A-B tilt (86.26 deg) equals the second
# principal angle printed in section 12, because for B the NAMED eigenplane is
# PC1-PC2 -- the same plane view A draws. The A-C tilt (88.64 deg) does NOT equal
# C's section-12 value (3.32 deg): that one is measured against span(V1,V3), a
# different plane. The two numbers are not in conflict; they answer different
# questions.
A_coords <- Xstd %*% V_A

add_other_plane <- function(V_other, fig_title, col = "#9467bd") {
  e1 <- as.numeric(t(V_A) %*% V_other[, 1])    # the other view's RC1, in A-coordinates
  e2 <- as.numeric(t(V_A) %*% V_other[, 2])    # ...................RC2
  n  <- as.numeric(t(V_A) %*% V_other[, 3])    # ... its residual axis = the plane's normal
  M  <- cbind(e1, e2, n)

  stopifnot(
    abs(sum(n^2) - 1) < 1e-8,                                   # n is a unit normal
    max(abs(crossprod(M) - diag(3))) < 1e-10,                   # orthonormal frame
    max(abs(A_coords %*% M - Xstd %*% V_other)) < 1e-10,        # A rotated by M IS the other view
    abs(sum(n * as.numeric(x_t %*% V_A))) < 1e-10               # plane passes through the target
  )

  # Size the quad to the cloud's spread along e1/e2, with a floor so it is always
  # visible even if the data barely extend in one of those directions.
  rad <- max(sqrt(rowSums(A_coords^2)))
  L1  <- max(max(abs(A_coords %*% e1)), 0.15 * rad) * 1.05
  L2  <- max(max(abs(A_coords %*% e2)), 0.15 * rad) * 1.05
  cor <- rbind(-L1 * e1 - L2 * e2, L1 * e1 - L2 * e2,
                L1 * e1 + L2 * e2, -L1 * e1 + L2 * e2)

  # TWO TILT ANGLES, NOT ONE. A plane through the origin in 3D has 2 degrees of
  # freedom (its unit normal, up to sign), so a single number cannot describe how it
  # sits. Writing the plane as  s3 = s1*s1_ax + s2*s2_ax  gives the two COMPONENT
  # tilts atan(s1) along PC1 and atan(s2) along PC2, and the dihedral angle is their
  # resultant:   tan^2(theta) = s1^2 + s2^2   (asserted below).
  #
  # The dihedral is the magnitude -- the minimum rotation that brings one plane onto
  # the other. The components say WHICH WAY it leans, which the magnitude hides. For
  # the target set in section 0, A-B is (-48.26 deg, +2.04 deg) -- almost a pure PC1
  # tilt, with the plane barely moving in the PC2 direction -- while A-C is
  # (+65.85 deg, +88.64 deg), leaning substantially both ways. Two planes, both
  # described by a single dihedral in the high 40s / high 80s, that sit quite
  # differently.
  #
  # Sign care: n is normalised to n[3] > 0 for REPORTING ONLY. The reconstruction
  # matrix M must keep the UNFLIPPED n, or A_coords %*% M no longer equals the other
  # view (a flipped copy makes that residual ~6 instead of ~1e-15).
  a_t   <- as.numeric(x_t %*% V_A)
  nrep  <- if (n[3] < 0) -n else n
  theta <- acos(min(1, abs(nrep[3])))
  if (abs(nrep[3]) < 1e-8) {
    s1 <- Inf; s2 <- Inf
  } else {
    s1 <- -nrep[1] / nrep[3]; s2 <- -nrep[2] / nrep[3]
    stopifnot(abs(tan(theta)^2 - (s1^2 + s2^2)) < 1e-8)
  }
  cat(sprintf(paste0("\n  %s\n",
      "    normal in A-coords       : (%.4f, %.4f, %.4f)\n",
      "    passes through target    : n . a_t = %.3e\n",
      "    tilt along PC1           : %+7.2f deg\n",
      "    tilt along PC2           : %+7.2f deg\n",
      "    dihedral (the resultant) : %7.2f deg   tan^2 = s1^2 + s2^2 = %.6f\n",
      "    A rotated by [e1 e2 n] reproduces the other view: max diff %.3e\n"),
      fig_title, n[1], n[2], n[3], sum(n * a_t),
      atan(s1) * 180 / pi, atan(s2) * 180 / pi, theta * 180 / pi, s1^2 + s2^2,
      max(abs(A_coords %*% M - Xstd %*% V_other))))

  # THE FREE FAMILY. "Contains the target" is ONE equation, s1*a1 + s2*a2 = a3, on a
  # plane that has two degrees of freedom -- so a whole 1-parameter family of tilts
  # is admissible and the Procrustes construction picks one of them arbitrarily
  # (section 15c). Sweeping that line shows how much better, or worse, the drawn
  # plane could have been for the same target.
  swp <- if (abs(a_t[2]) >= abs(a_t[1])) {
    g <- seq(-60, 60, length.out = 24001); cbind(g, (a_t[3] - g * a_t[1]) / a_t[2])
  } else {
    g <- seq(-60, 60, length.out = 24001); cbind((a_t[3] - g * a_t[2]) / a_t[1], g)
  }
  dih <- atan(sqrt(swp[, 1]^2 + swp[, 2]^2)) * 180 / pi
  cat(sprintf(paste0("    admissible dihedrals for this target: %.2f to %.2f deg",
                     "  (Procrustes gave %.2f)\n"),
              min(dih), max(dih), theta * 180 / pi))
  # Worth noticing across the two calls: for this target the pinned (1,2) plane comes
  # out at the SHALLOW end of that range and the selected pair's plane at the STEEP
  # end. How steeply a plane is tilted therefore says nothing about which pair
  # bl_find_local_cf() selects -- selection is on the Mahalanobis distance from the
  # target to its counterfactual (section 18e), and here the two run opposite ways.
  # A shallower plane is easier to read; it is not a better answer.

  # Built ON TOP of view A's 3D figure, so its points, green surface, black boundary,
  # loading arrows and target marker are all inherited unchanged.
  fig <- view_A$fig3d |>
    # 1) the plane itself: two triangles over the four corners (mesh3d rather than
    #    add_surface -- the plane is near-vertical here, so z = f(x,y) degenerates)
    add_trace(type = "mesh3d",
              x = cor[, 1], y = cor[, 2], z = cor[, 3],
              i = c(0, 0), j = c(1, 2), k = c(2, 3),
              facecolor = rep(col, 2), opacity = 0.25,
              hoverinfo = "skip", showlegend = FALSE) |>
    # 2) the in-plane axes -- exactly the other view's RC1 and RC2
    add_trace(x = c(-L1, L1) * e1[1], y = c(-L1, L1) * e1[2], z = c(-L1, L1) * e1[3],
              type = "scatter3d", mode = "lines+text",
              text = c("", "RC1"), textposition = "top center",
              textfont = list(size = 14, color = col),
              line = list(color = col, width = 6), showlegend = FALSE) |>
    add_trace(x = c(-L2, L2) * e2[1], y = c(-L2, L2) * e2[2], z = c(-L2, L2) * e2[3],
              type = "scatter3d", mode = "lines+text",
              text = c("", "RC2"), textposition = "top center",
              textfont = list(size = 14, color = col),
              line = list(color = col, width = 6), showlegend = FALSE) |>
    # 3) the normal -- "rotate until this points straight up"
    add_trace(x = c(0, n[1]) * rad * 0.6, y = c(0, n[2]) * rad * 0.6,
              z = c(0, n[3]) * rad * 0.6,
              type = "scatter3d", mode = "lines+text",
              text = c("", "normal (RC3)"), textposition = "top center",
              textfont = list(size = 13, color = col),
              line = list(color = col, width = 4, dash = "dot"), showlegend = FALSE)

  # 4) the other view's decision boundary. Its segments are (s1,s2) at s3 = 0, so in
  #    A-coordinates they are s1*e1 + s2*e2.
  list(fig = fig, e1 = e1, e2 = e2, n = n)
}

add_boundary_in_plane <- function(fig, cls, e1, e2) {
  for (cl in cls) {
    xy <- cbind(cl$x, cl$y) %*% rbind(e1, e2)      # -> A-coordinates
    fig <- fig |> add_trace(x = xy[, 1], y = xy[, 2], z = xy[, 3],
                            type = "scatter3d", mode = "lines",
                            line = list(width = 6, color = "red"),   # in-plane decision boundary
                            name = "other view's boundary", showlegend = FALSE)
  }
  fig
}

cat("\n--- 17) B and C shown inside view A ---")

ab <- add_other_plane(V_B, "A-B plane")
fig_AB_plane <- add_boundary_in_plane(ab$fig, view_B$cls, ab$e1, ab$e2) |>
  layout(title = "A-B plane - view A carrying B's biplot plane (rotated from (1,2)); tilt it horizontal to see B")

ac <- add_other_plane(V_C, sprintf("A-C plane (optimal pair (%d,%d))",
                                   best_pair[1], best_pair[2]))
fig_AC_plane <- add_boundary_in_plane(ac$fig, view_C$cls, ac$e1, ac$e2) |>
  layout(title = sprintf(paste0("A-C plane - view A carrying C's biplot plane ",
                                "(rotated from (%d,%d), optimal); tilt it horizontal to see C"),
                         best_pair[1], best_pair[2]))

# Run these in RStudio. In each, the violet plane should pass through the magenta
# target; rotate until the violet dotted normal points at you and the plane is
# face-on, and you are looking at view B (or C) -- with the black curve in the plane
# being that view's decision boundary.
fig_AB_plane
fig_AC_plane

# --- 18) THE 2D BIPLOT IS A SLICE, NOT A VIEW -------------------------------
# ----------------------------------------------------------------------------
# The question this answers: why do plot(loc_12) and plot(loc_op) look so
# different -- different local decision boundaries -- when the 3D figures look
# like the same object?
#
# Because they are different KINDS of object.
#
#   THE 3D FIGURE genuinely is the same object each time. With 3 variables the
#   basis is a lossless 3x3 rotation, so views A, B and C hold the same p ~ 0.5
#   cloud, the same 221 points and the same surface, re-expressed. Rotating only
#   moves the camera. (The green surface is re-fitted per frame, so it is a
#   re-regression of one point set rather than pixel-identical, but the object
#   underneath is one surface.)
#
#   THE 2D BIPLOT is a SLICE through that surface. plot(bl_local) draws where the
#   model crosses 0.5 WITHIN THE DRAWN PLANE, with the off-plane coordinate pinned
#   at 0. Two planes through the same target cut one surface along two genuinely
#   different curves. Neither is wrong; they are different cross-sections.
#
# And there is a geometric constraint that ties them together. Both planes contain
# the ORIGIN (a rotation has no translation) and the TARGET. Two distinct planes
# meet in exactly one line, so their intersection is FORCED to be span(x_t) -- the
# target's own ray. The surface crosses that ray at one point, and that point must
# therefore lie on BOTH 2D boundaries. The figure below shows exactly that: two
# curves on one surface, pinned together where the ray pierces it.
#
# The consequence is not cosmetic: a different slice gives a different nearest
# boundary point, hence a different counterfactual and different advice. That is
# why the pair search exists, and why the plane being under-determined
# (section 15c) matters.
if (CHECK_PKG) {

  to_feature <- function(lo) {          # a 2D biplot curve -> standardised feature space
    stopifnot(length(lo$ct_local) > 0L)
    M <- do.call(rbind, lapply(lo$ct_local, function(s) cbind(s$x, s$y)))
    list(xy = M, X = M %*% t(lo$Vr_rot))
  }
  s12 <- to_feature(loc_12)             # the PC1/2-pinned slice
  sop <- to_feature(loc_op)             # the optimal-pair slice

  # ---- 18a) both curves ARE p = 0.5 curves of the same model ---------------
  cat("\n--- 18a) Both 2D boundaries are genuine p = 0.5 curves ---\n")
  pfun <- function(X) { X <- as.data.frame(X); colnames(X) <- feat_names
                        predict(gam_fit, newdata = X, type = "response") }
  p12 <- pfun(s12$X); pop <- pfun(sop$X)
  cat(sprintf("  pinned (1,2) slice : %4d pts | p in [%.5f, %.5f] | max|p-0.5| = %.2e\n",
              nrow(s12$X), min(p12), max(p12), max(abs(p12 - 0.5))))
  cat(sprintf("  optimal (%d,%d) slice: %4d pts | p in [%.5f, %.5f] | max|p-0.5| = %.2e\n",
              pkg_pair[1], pkg_pair[2], nrow(sop$X), min(pop), max(pop), max(abs(pop - 0.5))))
  stopifnot(max(abs(p12 - 0.5)) < 5e-3, max(abs(pop - 0.5)) < 5e-3)
  cat("  Neither is wrong -- they are different cross-sections of one surface.\n")

  # ---- 18b) ... but they are different curves ------------------------------
  nn_dist <- function(A, B) apply(A, 1L, function(a) sqrt(min(colSums((t(B) - a)^2))))
  d_ab <- nn_dist(s12$X, sop$X); d_ba <- nn_dist(sop$X, s12$X)
  span <- max(dist(Xstd))
  cat("\n--- 18b) They are genuinely DIFFERENT curves ---\n")
  cat(sprintf(paste0("  nearest-neighbour distance between them: min %.4f | median %.4f |\n",
                     "  mean %.4f | max %.4f (Hausdorff %.4f), against a data span of %.3f\n"),
              min(d_ab), median(d_ab), mean(d_ab), max(d_ab),
              max(max(d_ab), max(d_ba)), span))
  cat(sprintf("  -> they separate by up to %.0f%% of the data's extent.\n",
              100 * max(d_ab) / span))

  # ---- 18c) but they MUST cross where the surface meets the target's ray ---
  cat("\n--- 18c) The one point they must share ---\n")
  nB <- as.numeric(t(V_A) %*% V_B[, 3]); nC <- as.numeric(t(V_A) %*% V_C[, 3])
  cx <- c(nB[2]*nC[3] - nB[3]*nC[2], nB[3]*nC[1] - nB[1]*nC[3], nB[1]*nC[2] - nB[2]*nC[1])
  cx <- cx / sqrt(sum(cx^2))
  ux <- x_t / sqrt(sum(x_t^2))
  cos_ray <- abs(sum(cx * as.numeric(t(V_A) %*% ux)))
  cat(sprintf("  planes intersect in the target's own ray? |cos| = %.12f\n", cos_ray))
  stopifnot(abs(cos_ray - 1) < 1e-8)

  # where does the surface cross that ray?
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
    # k0 is the scan point CLOSEST to 0.5, which need not be the FIRST crossing: the
    # surface can meet the ray more than once. n_cross, printed below, says how many
    # times it does for the target actually set.
    cat(sprintf("  surface crosses the ray %d time(s); closest scan point at t = %.5f -> (%+.4f, %+.4f, %+.4f), p = %.6f\n",
                n_cross, t_scan[k0], x_cross[1], x_cross[2], x_cross[3], p_ray[k0]))
    cat(sprintf("  (the target is at t = 1, so the crossing lies %s it)\n",
                if (t_scan[k0] < 1) "between the origin and" else "beyond"))
    d12 <- min(sqrt(colSums((t(s12$X) - x_cross)^2)))
    dop <- min(sqrt(colSums((t(sop$X) - x_cross)^2)))
    cat(sprintf("  distance from that point to each curve: pinned %.5f | optimal %.5f\n", d12, dop))
    cat(sprintf("  contour grid step, for scale: %.5f -> both pass through it\n",
                diff(loc_12$xseq)[1]))
    stopifnot(d12 < 3 * diff(loc_12$xseq)[1], dop < 3 * diff(loc_12$xseq)[1])
  }

  # ---- 18d) the in-plane axes differ too (a cosmetic second reason) --------
  cat("\n--- 18d) The panels are also drawn on different in-plane axes ---\n")
  cat(sprintf("  RC1 differs by %.2f deg | RC2 by %.2f deg\n",
      acos(min(1, abs(sum(loc_12$Vr_rot[,1] * loc_op$Vr_rot[,1])))) * 180/pi,
      acos(min(1, abs(sum(loc_12$Vr_rot[,2] * loc_op$Vr_rot[,2])))) * 180/pi))
  cat(sprintf("  target at (%.4f, %.4f) vs (%.4f, %.4f) -- norms %.4f and %.4f, both = ||x_t||\n",
              loc_12$Z_target[1], loc_12$Z_target[2],
              loc_op$Z_target[1], loc_op$Z_target[2],
              sqrt(sum(loc_12$Z_target^2)), sqrt(sum(loc_op$Z_target^2))))
  cat("  So even shared structure appears rotated between the two panels.\n")

  # ---- 18e) different slice => different counterfactual --------------------
  cat("\n--- 18e) Why this is not cosmetic ---\n")
  for (nm in c("pinned (1,2)", "optimal")) {
    lo <- if (nm == "pinned (1,2)") loc_12 else loc_op
    cat(sprintf("  %-13s B_x = (%s) | dist_z %.4f | sqrt(mah) %.4f\n", nm,
        paste(sprintf("%+.3f", as.numeric(lo$B_x)), collapse = ", "),
        lo$dist_z, sqrt(lo$dist_mahalanobis)))
  }
  cat(sprintf("  the two counterfactuals are %.4f apart in feature space.\n",
              sqrt(sum((as.numeric(loc_12$B_x) - as.numeric(loc_op$B_x))^2))))
  cat("  Different slice -> different nearest boundary point -> different advice.\n")

  # ---- 18f) the figure: BOTH planes and BOTH slices, on one surface --------
  # This is fig_AB_plane and fig_AC_plane merged into a single picture, so the two
  # planes can be compared directly instead of by flicking between figures.
  #   plane B (pinned (1,2)) : PURPLE     with its slice in RED
  #   plane C (optimal)      : LIGHT BLUE with its slice in DARK BLUE
  # Built on view A's 3D figure, so the points, green surface, loading arrows and
  # target marker come free. The RC1/RC2 axes and the plane normals are deliberately
  # left out -- section 17 carries those per plane, and six more lines here would
  # make the picture unreadable.
  toA <- function(sl, lo) sl$xy %*% t(lo$Vr_rot) %*% V_A   # curve -> A-coordinates
  A12 <- toA(s12, loc_12); Aop <- toA(sop, loc_op)
  A_coords <- Xstd %*% V_A
  rad <- max(sqrt(rowSums(A_coords^2)))
  ray <- as.numeric(t(V_A) %*% ux)

  # Same quad construction as add_other_plane() in section 17: size the plane to the
  # cloud's spread along its own in-plane directions, with a floor so it stays visible.
  quad <- function(V_other) {
    e1 <- as.numeric(t(V_A) %*% V_other[, 1])
    e2 <- as.numeric(t(V_A) %*% V_other[, 2])
    L1 <- max(max(abs(A_coords %*% e1)), 0.15 * rad) * 1.05
    L2 <- max(max(abs(A_coords %*% e2)), 0.15 * rad) * 1.05
    rbind(-L1 * e1 - L2 * e2, L1 * e1 - L2 * e2,
           L1 * e1 + L2 * e2, -L1 * e1 + L2 * e2)
  }
  cB <- quad(V_B); cC <- quad(V_C)

  add_quad <- function(fig, cor, col, nm) {
    add_trace(fig, type = "mesh3d",
              x = cor[, 1], y = cor[, 2], z = cor[, 3],
              i = c(0, 0), j = c(1, 2), k = c(2, 3),
              facecolor = rep(col, 2), opacity = 0.22,
              hoverinfo = "skip", name = nm, showlegend = TRUE)
  }

  # THE COUNTERFACTUAL ON EACH PLANE, and the measure that chose between them.
  #
  # bl_find_local_cf() works in two stages with two different metrics, and only the
  # second decides which pair is selected:
  #   WITHIN a plane -- the nearest boundary point to the target is found by plain
  #     EUCLIDEAN distance in the rotated 2D plane (dist_z).
  #   BETWEEN planes -- the winner is the pair with the smallest MAHALANOBIS distance
  #     from the target to that point, measured in feature space under the inverse
  #     pooled within-class covariance. It is stored SQUARED in dist_mahalanobis;
  #     the legend shows its square root, which is the distance itself.
  # The legend therefore carries the Mahalanobis value, because that is what "optimal"
  # means here. Section 18e prints both metrics side by side; they need not agree.
  cf_lab <- function(lo, nm) {
    if (is.null(lo$B_x)) return(NA_character_)
    if (is.na(lo$dist_mahalanobis))
      sprintf("CF | %s - Euclidean %.4f (Mahalanobis unavailable)", nm, lo$dist_z)
    else
      sprintf("CF | %s - Mahalanobis %.4f", nm, sqrt(lo$dist_mahalanobis))
  }
  sel_12 <- !is.na(loc_12$dist_mahalanobis) && !is.na(loc_op$dist_mahalanobis) &&
            loc_12$dist_mahalanobis < loc_op$dist_mahalanobis
  lab_12 <- paste0(cf_lab(loc_12, "pinned (1,2)"),
                   if (isTRUE(sel_12)) "  <- selected (smallest)" else "")
  lab_op <- paste0(cf_lab(loc_op, sprintf("optimal (%d,%d)", pkg_pair[1], pkg_pair[2])),
                   if (!isTRUE(sel_12)) "  <- selected (smallest)" else "")

  add_cf <- function(fig, lo, col, lab) {          # CF marker + target-to-CF connector
    if (is.null(lo$B_x)) return(fig)
    b <- as.numeric(as.numeric(lo$B_x) %*% V_A)    # counterfactual, in A-coordinates
    fig |>
      add_trace(x = c(a_tA[1], b[1]), y = c(a_tA[2], b[2]), z = c(a_tA[3], b[3]),
                type = "scatter3d", mode = "lines",
                line = list(color = col, width = 4, dash = "dot"),
                hoverinfo = "skip", showlegend = FALSE) |>
      add_markers(x = b[1], y = b[2], z = b[3],
                  marker = list(size = 9, color = col, symbol = "circle"),
                  name = lab, showlegend = TRUE)
  }
  a_tA <- as.numeric(x_t %*% V_A)                  # the target, in A-coordinates

  fig_slices <- view_A$fig3d |>
    add_quad(cB, "purple",    "plane: pinned (1,2)") |>
    add_quad(cC, "blue", sprintf("plane: optimal (%d,%d)", pkg_pair[1], pkg_pair[2])) |>
    add_trace(x = c(-rad, rad) * ray[1], y = c(-rad, rad) * ray[2],
              z = c(-rad, rad) * ray[3],
              type = "scatter3d", mode = "lines",
              line = list(color = "#666666", width = 4, dash = "dash"),
              name = "target's ray (where the planes meet)", showlegend = TRUE) |>
    add_trace(x = A12[, 1], y = A12[, 2], z = A12[, 3],
              type = "scatter3d", mode = "lines",
              line = list(color = "red", width = 6),
              name = "slice: pinned (1,2)", showlegend = TRUE) |>
    add_trace(x = Aop[, 1], y = Aop[, 2], z = Aop[, 3],
              type = "scatter3d", mode = "lines",
              line = list(color = "darkblue", width = 6),
              name = sprintf("slice: optimal (%d,%d)", pkg_pair[1], pkg_pair[2]),
              showlegend = TRUE)
  fig_slices <- add_cf(fig_slices, loc_12, "red",      lab_12)
  fig_slices <- add_cf(fig_slices, loc_op, "darkblue", lab_op)
  if (!is.null(x_cross)) {
    xcA <- as.numeric(t(V_A) %*% x_cross)
    fig_slices <- fig_slices |>
      add_markers(x = xcA[1], y = xcA[2], z = xcA[3],
                  marker = list(size = 10, color = "black", symbol = "circle"),
                  name = "both curves pass here", showlegend = TRUE)
  }
  fig_slices <- fig_slices |>
    layout(title = paste("Both rotated planes on ONE surface - each 2D biplot is a",
                         "cross-section, not a view.<br>",
                         "<sup>Pair selected by the smallest Mahalanobis distance from",
                         "the target to its counterfactual.</sup>"),
           showlegend = TRUE)

  # Run in RStudio: two translucent planes cutting one green surface, each carrying
  # its own boundary curve, and the two curves meeting where the dashed ray -- the
  # line the planes share -- pierces the surface.
  fig_slices
}
