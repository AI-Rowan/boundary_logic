# Originated in the `Explain XAI code` project's sandbox as
# `R/In-development/17_1a_toy-2d-to-1d-rotation.R` (snapshot 2026-09-10). This repo
# is now this script's home; edit it here.

# ============================================================================
# The Procrustes rotation at minimum scale — 2 variables reduced to a 1D biplot
# ----------------------------------------------------------------------------
# Supersedes "A3 Rotated Biplot.R", from which the data and the panel layout are
# taken. Redone in line with script 17_1: RC naming, trace-based RSS, and the
# adequacy / predictivity block.
#
# Five points in two variables are reduced to a ONE-dimensional biplot (r = 1),
# and the basis is rotated so that target point 5 lies exactly on that line.
# Because there is only one component left to draw, this is the smallest setting
# in which the rotation can be shown at all -- which is why it is a good place to
# build intuition, and why that intuition needs the warning at the bottom.
#
# THE QUESTION THIS SCRIPT ANSWERS: rotate onto PC1, or onto PC2?
#
# The answer is that IT MAKES NO DIFFERENCE. The two options produce the same
# rotation matrix A, the same retained axis, the same adequacy, the same RSS and
# the same sample predictivity. Only rho differs (1.134 vs 2.120), and rho is
# computed and then deliberately discarded -- by the thesis and by the package
# alike -- so nothing that is drawn depends on it.
#
# WHY THEY COINCIDE. The padding cbind(Y %*% Vr, 0) puts the retained score in
# SLOT 1 whichever component is nominated, so the target's reduced representation
# points along e1 either way. The Procrustes cross-product is C = 2*u*w' with
# w = (u_keep, 0): the same rank-1 column space (u) and row space (e1) in both
# cases, differing only by the scalar u_keep. A positive scalar does not change
# an SVD's singular vectors, so A comes out identical.
#
# The geometric version of the same statement: the retained axis is the TARGET'S
# OWN DIRECTION (verified below for all five points). In 2D there is exactly one
# line through the origin that contains the target, so there is nothing to choose.
#
# WHAT "DISAPPEARS", PRECISELY. Not PC2. BOTH original components disappear: the
# retained axis is RC1 = the target's direction, which is neither PC1 nor PC2, and
# the dropped axis is RC2 = its perpendicular. Per-variable adequacy moves from
# (0.3205, 0.6795) to (0.0122, 0.9878) -- variable 1's representation nearly
# vanishes while variable 2's becomes near-total. This is the r = 1 form of the
# point 17_1 makes at r = 2: after the rotation the original PC labels stop
# applying, and the axes are rotated components.
#
# WHERE THIS INTUITION BREAKS -- read this before generalising from the picture.
#
#   setting            subspace DOF   minus "contains the target"   free
#   p=2, r=1 (here)          1                  -1                   0   unique
#   p=3, r=2 (17_1)          2                  -1                   1
#   p=4, r=2 (17_2)          4                  -1                   3
#
# This toy is the ONLY case in the family where the rotation is fully determined.
# From p = 3 on there is leftover freedom, the Procrustes SVD resolves it
# arbitrarily, and the resulting plane is not reproducible -- see 17_1 section 15c
# and biplot_rotation_and_pair_selection.md section 4.1.
#
# A residue of that freedom is visible even here: for points 1 and 4 the two
# options give DIFFERENT A matrices, yet the same retained axis. At p = 2 the
# degenerate subspace is one-dimensional, so the only freedom is the orientation
# of the DISCARDED axis -- which a 1D biplot never shows.
#
# CAVEAT ON THE THESIS NUMBERS. theta = 28.14 deg and rho = 1.134 reproduce the
# thesis's "28 degree rotation, rho = 1/1.13" (note section 2.5) exactly. But the
# thesis reports sample predictivity 0.72 -> 1 and RSS 0.28 -> 4 of 14, where this
# data gives 0.7776 -> 1 and 4.307 -> 4.927 of 11.4. These five points are
# therefore probably NOT the thesis's dataset; do not read the agreement in theta
# and rho as reproduction of the worked example.
#
# NOTE ON GRAPHICS. Base R, deliberately: a 5-point 2D example reads better as
# annotated static panels. Base-R plots do NOT render through Quarto's `#| file:`
# (CLAUDE.md section 6), so this would need a static export or a plotly rewrite
# before any article use.
# ============================================================================

rm(list = ls())

# --- 0) Data (A3's five points), centred but NOT scaled ---------------------
# Centring only, as A3 does. Consequently tr(X'X) = (n-1)*p does NOT apply here
# (that identity needs standardised data); total SS is 11.4, not 8. The other
# trace identity does hold and is asserted below: unrotated RSS == d[2]^2.
X.start <- matrix(c(1.5, 1.5,
                    4,   2,
                    2,   4,
                    4,   4,
                    3,   4), ncol = 2, byrow = TRUE)
X.st <- scale(X.start, center = TRUE, scale = FALSE)
p    <- 2L
i    <- 5L                       # target observation
n    <- nrow(X.st)

svdx <- svd(X.st); V <- svdx$v; D <- svdx$d
x    <- as.numeric(X.st[i, ])
u    <- as.numeric(x %*% V)      # the target in PC coordinates
stopifnot(sqrt(sum(x^2)) > 1e-10, min(D) > 1e-8)

XV2 <- X.st %*% V                # scores on both PCs
d2r <- function(r) r * 180 / pi

# --- 1) The rotation, with the retained component as a parameter ------------
# A3's construction, generalised: `keep` nominates which component becomes the
# 1D biplot. Everything else is A3 line for line.
a3_rot <- function(keep) {
  stopifnot(abs(u[keep]) > 1e-10)   # p=2 form of 17_1's ||w|| -> 0 warning:
                                    # a zero score makes C the zero matrix and A arbitrary
  X3 <- rbind(-X.st[i, , drop = FALSE], matrix(0, 1, p), X.st[i, , drop = FALSE])
  Y  <- X3 %*% V                                      # full-PC representation
  Z  <- X3 %*% V[, keep, drop = FALSE]                # reduced representation
  Z  <- cbind(Z, matrix(0, ncol = p - 1L, nrow = 3))  # padded back out to p columns
  C  <- t(Y) %*% Z
  s  <- svd(C)
  A  <- s$v %*% t(s$u)
  list(A    = A,
       Vrho = V %*% t(A),
       rho  = sum(diag(Z %*% A %*% t(Y))) / sum(diag(t(Z) %*% Z)),
       dsv  = s$d)
}

r1 <- a3_rot(1L)     # option 1: keep PC1
r2 <- a3_rot(2L)     # option 2: keep PC2

cat("=== The target ===\n")
cat(sprintf("  point %d, PC coords u = (%+.4f, %+.4f), ||u|| = %.4f\n",
            i, u[1], u[2], sqrt(sum(u^2))))
cat(sprintf("  angle to PC1 %.2f deg | angle to PC2 %.2f deg | sum %.2f deg (forced in 2D)\n",
            d2r(acos(abs(u[1]) / sqrt(sum(u^2)))),
            d2r(acos(abs(u[2]) / sqrt(sum(u^2)))),
            d2r(acos(abs(u[1]) / sqrt(sum(u^2)))) +
            d2r(acos(abs(u[2]) / sqrt(sum(u^2))))))

cat("\n=== Option 1 (keep PC1) vs option 2 (keep PC2) ===\n")
cat(sprintf("  Procrustes C singular values : keep PC1 %s | keep PC2 %s  -> rank 1 of %d\n",
            paste(format(r1$dsv, digits = 3), collapse = ", "),
            paste(format(r2$dsv, digits = 3), collapse = ", "), p))
cat(sprintf("  applied rotation acos(A[1,1]): keep PC1 %.2f deg | keep PC2 %.2f deg\n",
            d2r(acos(min(1, abs(r1$A[1, 1])))), d2r(acos(min(1, abs(r2$A[1, 1]))))))
cat(sprintf("  retained axis RC1            : keep PC1 (%.6f, %.6f) | keep PC2 (%.6f, %.6f)\n",
            r1$Vrho[1, 1], r1$Vrho[2, 1], r2$Vrho[1, 1], r2$Vrho[2, 1]))
cat(sprintf("  |cos| between the two axes   : %.10f\n",
            abs(sum(r1$Vrho[, 1] * r2$Vrho[, 1]))))
cat(sprintf("  rho (computed, NOT applied)  : keep PC1 %.6f | keep PC2 %.6f  <- the only difference\n",
            r1$rho, r2$rho))

cat("\n=== The retained axis is the target's own direction ===\n")
ux <- x / sqrt(sum(x^2))
cat(sprintf("  x/||x|| = (%.6f, %.6f)\n", ux[1], ux[2]))
for (k in 1:n) {
  xk <- as.numeric(X.st[k, ]); uk <- as.numeric(xk %*% V)
  X3k <- rbind(-X.st[k, , drop = FALSE], matrix(0, 1, p), X.st[k, , drop = FALSE])
  Yk  <- X3k %*% V
  Zk  <- cbind(X3k %*% V[, 1L, drop = FALSE], matrix(0, ncol = p - 1L, nrow = 3))
  sk  <- svd(t(Yk) %*% Zk); Vk <- V %*% t(sk$v %*% t(sk$u))
  cat(sprintf("  point %d: |cos(RC1, x)| = %.12f\n", k,
              abs(sum(Vk[, 1] * xk / sqrt(sum(xk^2))))))
}

# --- 2) Diagnostics: what the rotation costs --------------------------------
# 17_1's section-13 block, with J = diag(1, 0) because the biplot here is 1D.
# RSS uses the TRACE (A3's original summed every matrix entry, which is not a sum
# of squares and made the rotation look like it improved global fit).
J      <- diag(c(1, 0))
S      <- D^2
XtX    <- crossprod(X.st)
tot_SS <- sum(diag(XtX))

rc_diagnostics <- function(Vv, label) {
  # Uv uses the original D; J is diagonal so it commutes with D and the factors
  # cancel, leaving ||proj of x_k onto the kept axis||^2 / ||x_k||^2.
  Uv  <- X.st %*% Vv %*% solve(diag(D))
  adq <- diag(Vv %*% J %*% t(Vv)) / diag(Vv %*% t(Vv))
  axp <- diag(Vv %*% diag(S) %*% J %*% t(Vv)) / diag(Vv %*% diag(S) %*% t(Vv))
  smp <- diag(Uv %*% diag(S) %*% J %*% t(Uv)) / diag(Uv %*% diag(S) %*% t(Uv))
  rss <- tot_SS - sum(diag(XtX %*% Vv %*% J %*% t(Vv)))
  xh  <- max(abs(x %*% Vv[, 1, drop = FALSE] %*% t(Vv[, 1, drop = FALSE]) - x))
  cat(sprintf(paste0("  %-12s adequacy %s | axis-pred %s | RSS %7.4f of %.4f",
                     " | target SP %.6f | others %.4f | max|xhat-x| %.2e\n"),
              label, paste(sprintf("%.4f", adq), collapse = " "),
              paste(sprintf("%.4f", axp), collapse = " "),
              rss, tot_SS, smp[i], mean(smp[-i]), xh))
  list(adq = adq, axp = axp, smp = smp, rss = rss, xh = xh)
}

cat("\n=== What the rotation costs (J = diag(1,0): the 1D biplot) ===\n")
dgA <- rc_diagnostics(V,       "unrotated")
dg1 <- rc_diagnostics(r1$Vrho, "keep PC1")
dg2 <- rc_diagnostics(r2$Vrho, "keep PC2")

cat(sprintf(paste0("\n  THE TRADE: target %d becomes EXACT (%.4f -> %.6f) while the other %d\n",
                   "  points get slightly worse (%.4f -> %.4f) and the residual rises\n",
                   "  (%.4f -> %.4f of %.4f). Variable %d's representation is what collapses\n",
                   "  (adequacy %.4f -> %.4f): the axes are no longer PC1 and PC2.\n"),
            i, dgA$smp[i], dg1$smp[i], n - 1L,
            mean(dgA$smp[-i]), mean(dg1$smp[-i]),
            dgA$rss, dg1$rss, tot_SS,
            which.min(dg1$adq), dgA$adq[which.min(dg1$adq)], min(dg1$adq)))

# --- 3) The equivalence, asserted -------------------------------------------
stopifnot(
  max(abs(r1$A - r2$A)) < 1e-12,                                   # same rotation
  abs(abs(sum(r1$Vrho[, 1] * r2$Vrho[, 1])) - 1) < 1e-12,          # same retained axis
  abs(abs(sum(r1$Vrho[, 1] * ux)) - 1) < 1e-12,                    # = the target's direction
  abs(dg1$rss - dg2$rss) < 1e-10,                                  # same cost
  max(abs(dg1$adq - dg2$adq)) < 1e-10,
  abs(dg1$smp[i] - 1) < 1e-8, abs(dg2$smp[i] - 1) < 1e-8,          # target becomes exact
  dg1$xh < 1e-10, dg2$xh < 1e-10, dgA$xh > 1e-3,                   # xhat = x only after
  abs(dgA$rss - D[2]^2) < 1e-10                                    # unrotated RSS is d2^2
)
cat("\n  All equivalence assertions pass: the two options are the same biplot.\n")

# --- 4) Degrees of freedom: why this toy is the only unambiguous case -------
cat("\n=== Degrees of freedom ===\n")
# CORRECTED 2026-09-10. This table previously ran subspace / constraint / free and
# gave "free = 3" at p = 4, citing the note's "3-dimensional manifold". That 3 is
# the count of free ROTATIONS, not of free PLANES, and the two are different
# numbers: containing a fixed line is TWO conditions on Gr(2,4), not one, so the
# plane family has dimension 2. Measured numerically in 17_3 section 4 (tangent
# rank 2) and in 17_2 section 8a. The difference between the columns is the
# leftover freedom in the frame once the plane is fixed.
dof <- data.frame(
  setting     = c("p=2, r=1  (this script)", "p=3, r=2  (17_1)", "p=4, r=2  (17_2)"),
  planes      = c(0L, 1L, 2L),
  rotations   = c(1L, 2L, 3L),
  consequence = c("plane unique -- both options coincide",
                  "1 free plane direction; plane not reproducible",
                  "2 free plane directions; 3-dim manifold of rotations")
)
print(dof, row.names = FALSE)
cat("  'planes' counts admissible 2-planes (or lines, at r = 1) through the target;\n")
cat("  'rotations' counts the orthogonal matrices the Procrustes step may return.\n")
cat("  The rotation is fully determined ONLY here. Do not carry this picture's\n")
cat("  tidiness into p >= 3 -- see 17_1 section 15c and the note's section 4.1.\n")

# --- 5) Panels ---------------------------------------------------------------
op <- par(no.readonly = TRUE)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
cols <- c(1, 1, 1, 1, 2)          # target (point 5) in red
lim  <- c(-2.6, 2.6)

## A) the data in variable space
plot(X.st, asp = 1, col = cols, pch = 16, cex = 1.5, xlim = lim, ylim = lim,
     xlab = "Variable 1", ylab = "Variable 2", main = "A) Two-variable data")
text(X.st[, 1], X.st[, 2], labels = seq_len(n), pos = 3, cex = 0.9)
abline(h = 0, v = 0, col = "grey90")

## B) unrotated biplot, keeping PC1 -- the first option's starting point
plot(XV2, asp = 1, col = cols, pch = 16, cex = 1.5, xlim = lim, ylim = lim,
     xlab = "PC1", ylab = "PC2", main = "B) 1D biplot on PC1 (unrotated)")
for (k in 1:p) {
  arrows(0, 0, V[k, 1], V[k, 2], col = "seagreen", length = 0.1)
  text(V[k, 1], V[k, 2], labels = paste0("V", k), pos = 2, cex = 0.9)
}
abline(h = 0, lwd = 2, col = "dodgerblue")
segments(XV2[, 1], 0, XV2[, 1], XV2[, 2], col = "grey70")
points(XV2[, 1], rep(0, n), pch = 15, cex = 1.4, col = cols)
text(XV2[, 1], -0.25, labels = seq_len(n), cex = 0.8, col = "grey40")

## C) unrotated biplot, keeping PC2 -- the second option's starting point
plot(XV2, asp = 1, col = cols, pch = 16, cex = 1.5, xlim = lim, ylim = lim,
     xlab = "PC1", ylab = "PC2", main = "C) 1D biplot on PC2 (unrotated)")
for (k in 1:p) {
  arrows(0, 0, V[k, 1], V[k, 2], col = "seagreen", length = 0.1)
  text(V[k, 1], V[k, 2], labels = paste0("V", k), pos = 2, cex = 0.9)
}
abline(v = 0, lwd = 2, col = "dodgerblue")
segments(0, XV2[, 2], XV2[, 1], XV2[, 2], col = "grey70")
points(rep(0, n), XV2[, 2], pch = 15, cex = 1.4, col = cols)
text(0.25, XV2[, 2], labels = seq_len(n), cex = 0.8, col = "grey40")

## D) the rotation that is actually applied
theta <- acos(min(1, abs(r1$A[1, 1])))
plot(XV2, asp = 1, col = cols, pch = 16, cex = 1.5, xlim = lim, ylim = lim,
     xlab = "PC1", ylab = "PC2", main = "D) The rotation applied")
abline(h = 0, col = "grey85")
for (k in 1:p) {
  arrows(0, 0, V[k, 1], V[k, 2], col = "seagreen", length = 0.1)
  text(V[k, 1], V[k, 2], labels = paste0("V", k), pos = 2, cex = 0.9)
}
arrows(0, 0, XV2[i, 1], XV2[i, 2], length = 0.06)          # target's vector
arc <- seq(0, theta, length.out = 200)
lines(cos(arc) * sqrt(sum(u^2)), sin(arc) * sqrt(sum(u^2)), col = "grey40")
points(sqrt(sum(u^2)), 0, pch = 15, cex = 1.4, col = 2)    # where it lands (full length)
points(XV2[i, 1], 0, pch = 15, cex = 1.2, col = "grey50")  # orthogonal projection (shorter)
text(sqrt(sum(u^2)) / 2, 0.18,
     labels = sprintf("theta = %.1f deg", d2r(theta)), cex = 0.9)
text(-1.9, -2.2, labels = sprintf("rho = %.3f (not applied)", r1$rho),
     cex = 0.85, col = "grey30", pos = 4)

## E) the rotated component (RC) biplot
XVrho <- X.st %*% r1$Vrho
plot(XV2, asp = 1, col = cols, pch = 16, cex = 1.5, xlim = lim, ylim = lim,
     xlab = "PC1", ylab = "PC2", main = "E) RC biplot (keep PC1)")
abline(h = 0, col = "grey90")
L  <- max(abs(XVrho[, 1])) * 1.25
uu <- r1$Vrho[, 1] / sqrt(sum(r1$Vrho[, 1]^2))
# the retained axis, drawn in PC coordinates
ax <- as.numeric(t(V) %*% uu)
segments(-L * ax[1], -L * ax[2], L * ax[1], L * ax[2], lwd = 2, col = "dodgerblue")
s  <- as.numeric(XV2 %*% ax)
segments(XV2[, 1], XV2[, 2], s * ax[1], s * ax[2], col = "grey70")
points(s * ax[1], s * ax[2], pch = 15, cex = 1.4, col = cols)
for (k in 1:p) {
  arrows(0, 0, V[k, 1], V[k, 2], col = "seagreen", length = 0.1)
  text(V[k, 1], V[k, 2], labels = paste0("V", k), pos = 2, cex = 0.9)
}
text(XV2[, 1], XV2[, 2], labels = seq_len(n), pos = 1, cex = 0.8, col = "grey20")
text(-2.4, 2.3, labels = "target lies ON the line", cex = 0.85, col = 2, pos = 4)

## F) the punchline -- both options give the same line
plot(XV2, asp = 1, col = cols, pch = 16, cex = 1.5, xlim = lim, ylim = lim,
     xlab = "PC1", ylab = "PC2", main = "F) Both options coincide")
abline(h = 0, col = "grey90")
ax1 <- as.numeric(t(V) %*% (r1$Vrho[, 1] / sqrt(sum(r1$Vrho[, 1]^2))))
ax2 <- as.numeric(t(V) %*% (r2$Vrho[, 1] / sqrt(sum(r2$Vrho[, 1]^2))))
segments(-L * ax1[1], -L * ax1[2], L * ax1[1], L * ax1[2],
         lwd = 5, col = "dodgerblue")                       # from keeping PC1
segments(-L * ax2[1], -L * ax2[2], L * ax2[1], L * ax2[2],
         lwd = 2, col = "orange", lty = 2)                  # from keeping PC2
points(XV2[i, 1], XV2[i, 2], pch = 16, cex = 1.8, col = 2)
legend("topleft", bty = "n", cex = 0.75,
       legend = c("rotate onto PC1", "rotate onto PC2 (dashed)"),
       col = c("dodgerblue", "orange"), lty = c(1, 2), lwd = c(5, 2))
text(-2.4, -1.9, cex = 0.8, col = "grey20", pos = 4,
     labels = "1 DOF - 1 constraint = 0:\nthe line through origin and\ntarget is unique")

par(op)

cat("\nPanels drawn. Panel F is the demonstration: the two lines superimpose.\n")
