# Originated in the `Explain XAI code` project's sandbox as
# `R/In-development/17_3_minimal-rotation-verification.R` (snapshot 2026-09-10).
# This repo is now this script's home; edit it here.

# ============================================================================
# Verification harness for the reference note's §6 recommendation
# ----------------------------------------------------------------------------
# WHY THIS SCRIPT EXISTS
#
# §6 of biplot_rotation_and_pair_selection.md
# recommends replacing the rank-1 Procrustes rotation with a closed-form minimal
# (Givens) rotation, and replacing the six-eigenvector-pair sample with a search
# over the family of planes containing the target. Every number §6 quotes was
# measured -- and §4.1 declined the change pending exactly that measurement.
#
# This script IS that measurement. It is the harness the package project needs in
# order to re-run the acceptance tests before adopting anything, and it is why
# §6 is a proposal with evidence rather than an assertion.
#
# WHAT IS ASSERTED AND WHAT IS ONLY REPORTED
#
#   ASSERTED -- the geometry. Exactness, the (0, theta) principal angles, bitwise
#     determinism, the plane dimension, the monotonicity of the family search, and
#     the thesis reproduction. These are properties of the construction and must
#     hold on any machine. Each carries its expected value in a comment, so a
#     future drift is diagnosable rather than just a failure.
#   REPORTED ONLY -- anything timed, and the within-pair cost. Wall clock is
#     machine-dependent; asserting it would make this script fail on a different
#     laptop for no methodological reason.
#
# SECTIONS
#   1  min_rotation()          the closed-form construction (§6.2)
#   2  exactness + labels      the headline table (§6.2)
#   3  determinism             (§6.5 test 2)
#   4  plane dimension         the correction to §5.1's table (§6.1)
#   5  within-pair cost        (§6.2's supporting measurement)
#   6  ray_scan()              the in-plane objective (§6.3)
#   7  p = 3 exhaustive sweep  (§6.3's table)
#   8  p = 4 family lattice    (§6.3)
#   9  cost                    (§6.3's cost paragraph)
#  10  thesis reproduction     (§6.2)
#
# Console output only -- no figures. Runs in roughly a minute.
# ============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(boundarylogic)
  library(mgcv)
})

prin_ang <- function(Q1, Q2)
  acos(pmin(1, pmax(-1, svd(t(Q1) %*% Q2)$d))) * 180 / pi
Pj <- function(M) M %*% t(M)          # plane projector

gen_pairs <- function(p) {            # the package's colexicographic order
  out <- list()
  for (j in 2L:p) for (ii in 1L:(j - 1L)) out[[length(out) + 1L]] <- c(ii, j)
  do.call(rbind, out)
}

# ---- 0) The two data setups ------------------------------------------------
# p = 4 row 49 (the note's §4.1 target) and p = 3 row 6 (§5's workhorse target).
# Mirrors 17_2 §1-5 and 17_1 §1-5: standardise in place so the package's "raw"
# feature space IS this script's Xstd space, then wrap the GAM.
build_setup <- function(feats, target_row, k_smooth = 8L, ttl = "verify") {
  loan_raw <- read.csv(
    system.file("extdata", "loan_data.csv", package = "boundarylogic")
  )
  lf <- loan_raw[loan_raw$previous_loan_defaults_on_file == "No", ]
  lf <- lf[, c(feats, "loan_status")]
  lf$credit_score <- as.numeric(lf$credit_score)

  bl_dat <- bl_prepare_data(data = lf, class_col = "loan_status",
                            feature_cols = feats, train_fraction = 0.01,
                            seed = 121L, hull_fraction = 1)
  dat <- bl_dat$train_data
  sc  <- scale(dat[, feats], center = TRUE, scale = TRUE)
  Xc  <- attr(sc, "scaled:center"); Xs <- attr(sc, "scaled:scale")
  Xstd <- sc

  dg <- data.frame(Xstd, class = dat$class); names(dg)[seq_along(feats)] <- feats
  fml <- as.formula(paste("class ~",
    paste(sprintf("s(%s, k = %d)", feats, k_smooth), collapse = " + ")))
  gam_fit <- gam(fml, data = dg, family = binomial())

  bl_dat$train_data[, feats] <- Xstd
  bl_dat$test_data[,  feats] <- scale(bl_dat$test_data[, feats],
                                      center = Xc, scale = Xs)
  bl_dat <- bl_set_scaling(bl_dat, center = Xc, scale = Xs, method = "z-score")
  bl_mod <- bl_wrap_model(gam_fit, "GAM", bl_dat$var_names,
                          train_data = bl_dat$train_data)
  bl_res <- bl_build_result(bl_data = bl_dat, bl_model = bl_mod,
                            method = "PCA", title = ttl)

  pfun <- function(X) {
    X <- as.data.frame(X); colnames(X) <- feats
    predict(gam_fit, newdata = X, type = "response")
  }
  list(feats = feats, Xstd = Xstd, p = ncol(Xstd), i = target_row,
       x_t = as.numeric(Xstd[target_row, ]), V = bl_res$V,
       cutoff = bl_res$cutoff, Minv = bl_res$metric_inv,
       bl_res = bl_res, pfun = pfun,
       tgt = bl_select_target(bl_res, target = bl_res$train_data[target_row, ]))
}

cat("=== 0) Building the two setups ===\n")
S4 <- build_setup(c("loan_percent_income", "credit_score",
                    "loan_int_rate", "loan_amnt"), 49L, ttl = "p4")
S3 <- build_setup(c("loan_percent_income", "credit_score",
                    "loan_int_rate"),          6L,  ttl = "p3")
cat(sprintf("  p = 4, target row %d | p = 3, target row %d\n", S4$i, S3$i))

# The package's own answer, taken NOW -- as the FIRST bl_find_local_cf call on each
# setup, before anything else in this script has touched it. That ordering matters:
# section 8b measures a full search returning a DIFFERENT distance when it follows a
# max_pairs = 1L call on the same objects. Capturing the reference here keeps
# sections 7 and 8 comparing against the unperturbed figure.
pkg_ref4 <- bl_find_local_cf(S4$bl_res, S4$tgt, verbose = FALSE)$dist_mahalanobis
pkg_ref3 <- bl_find_local_cf(S3$bl_res, S3$tgt, verbose = FALSE)$dist_mahalanobis
cat(sprintf("  package reference (first calls): p = 4 %.6f | p = 3 %.6f\n",
            pkg_ref4, pkg_ref3))          # expect 1.147879 and 1.312277

# ---- 1) The proposed construction (§6.2) -----------------------------------
# Minimal Givens rotation carrying the named eigenplane onto a plane containing
# the target. R is a rotation by theta inside span(u, w) and the identity on its
# complement, so it is the SMALLEST rotation achieving exactness.
min_rotation <- function(x_t, V, pair) {
  p   <- length(x_t)
  Vr0 <- V[, pair, drop = FALSE]
  xP  <- as.numeric(Vr0 %*% (t(Vr0) %*% x_t))     # the target's shadow
  nx  <- sqrt(sum(x_t^2)); nP <- sqrt(sum(xP^2))

  # Edge case: the target is at the centroid -- every plane contains it and no
  # rotation is defined or needed.
  if (nx < 1e-12) {
    warning("target is at the origin; returning the unrotated eigenplane.",
            call. = FALSE)
    return(list(Vr = Vr0, R = diag(p), theta = 0))
  }
  # Edge case: the target is orthogonal to the eigenplane, so its shadow is zero
  # and u is undefined. Convention: rotate the FIRST named axis onto the target.
  if (nP < 1e-12) { u <- Vr0[, 1]; r <- x_t } else { u <- xP / nP; r <- x_t - xP }

  nr <- sqrt(sum(r^2))
  # Edge case: the target already lies in the eigenplane -- theta = 0, R = I.
  if (nr < 1e-12) return(list(Vr = Vr0, R = diag(p), theta = 0))

  w  <- r / nr
  ct <- min(1, nP / nx); st <- sqrt(max(0, 1 - ct^2))
  R  <- diag(p) + (ct - 1) * (outer(u, u) + outer(w, w)) +
        st * (outer(w, u) - outer(u, w))
  list(Vr = R %*% Vr0, R = R, theta = acos(ct) * 180 / pi)
}

# The CURRENT construction, for contrast: rank-1 Procrustes on the three points.
procrustes_plane <- function(x_t, V, pair) {
  p  <- length(x_t)
  Y  <- rbind(-x_t, rep(0, p), x_t)
  sv <- svd(t(Y %*% V) %*% cbind(Y %*% V[, pair], matrix(0, 3L, p - 2L)))
  (V %*% t(sv$v %*% t(sv$u)))[, 1:2, drop = FALSE]
}

# ---- 2) Exactness and label integrity, p = 4 (§6.2's table) ----------------
# The sharp result. Under the minimal rotation the plane SHARES a direction with
# the eigenplane its label names (first principal angle 0) and its only departure
# is the theta it was forced to move. Under the current construction it shares
# nothing -- which is the mechanism behind §5.2's overlap figures.
cat("\n=== 2) Exactness and label integrity (p = 4, row 49) ===\n")
prs4 <- gen_pairs(S4$p)
pn4  <- apply(prs4, 1L, function(r) sprintf("(%d,%d)", r[1], r[2]))
cat(sprintf("%8s %10s %11s %11s %12s %13s\n",
            "pair", "theta", "ang1_min", "ang2_min", "exactness", "ang1_current"))
for (k in seq_len(nrow(prs4))) {
  mr <- min_rotation(S4$x_t, S4$V, prs4[k, ])
  a  <- prin_ang(mr$Vr, S4$V[, prs4[k, ], drop = FALSE])
  ex <- sqrt(sum((S4$x_t - as.numeric(mr$Vr %*% (t(mr$Vr) %*% S4$x_t)))^2))
  ap <- prin_ang(procrustes_plane(S4$x_t, S4$V, prs4[k, ]),
                 S4$V[, prs4[k, ], drop = FALSE])
  cat(sprintf("%8s %10.4f %11.2e %11.4f %12.2e %13.4f\n",
              pn4[k], mr$theta, a[1], a[2], ex, ap[1]))
  # Expected: ang1 = 0.00e+00, ang2 = theta, exactness 1.1e-15 to 4.2e-15.
  # Expected ang1_current: 1.1111 / 5.1471 / 55.6173 / 17.2008 / 5.4437 / 51.3428
  stopifnot(ex < 1e-12, a[1] < 1e-8, abs(a[2] - mr$theta) < 1e-8)
}
cat("  ASSERTED: ang1 = 0 and ang2 = theta exactly, every pair.\n")
cat("  The rightmost column is the CURRENT construction's first angle -- nonzero,\n")
cat("  so its plane shares no direction with the eigenplane it is named after.\n")

# ---- 3) Determinism (§6.5 test 2) ------------------------------------------
cat("\n=== 3) Determinism ===\n")
h <- replicate(5L, paste(sprintf("%.17g",
       as.numeric(min_rotation(S4$x_t, S4$V, c(2L, 3L))$Vr)), collapse = "|"))
cat(sprintf("  basis bit-identical across 5 calls: %s\n", length(unique(h)) == 1L))
stopifnot(length(unique(h)) == 1L)
cat("  (The current construction fails the equivalent check: §5.3 records\n")
cat("   best_pair flipping between (1,2) and (1,3) across five identical calls.)\n")

# ---- 4) The plane dimension (§6.1) -----------------------------------------
# §5.1's table reports "free = 3" at p = 4, citing §4.1's "3-dimensional manifold
# of rotations". That 3 is right for ROTATIONS and wrong for PLANES: containing a
# fixed line is TWO conditions on Gr(2,4), not one, leaving 4 - 2 = 2. Measured
# here as the tangent rank of the planes-through-the-target family, rather than
# asserted algebraically.
cat("\n=== 4) Dimension of the planes-through-the-target family (p = 4) ===\n")
{
  p <- S4$p
  ray <- S4$x_t / sqrt(sum(S4$x_t^2))
  vr0 <- min_rotation(S4$x_t, S4$V, c(1L, 2L))$Vr
  d0  <- vr0[, 2] - vr0[, 1] * sum(vr0[, 1] * vr0[, 2])
  d0  <- d0 - ray * sum(d0 * ray); d0 <- d0 / sqrt(sum(d0^2))
  Q   <- cbind(ray, d0)
  cmp <- svd(diag(p) - Pj(Q))$u[, 1:(p - 2), drop = FALSE]
  eps <- 1e-6; P0 <- Pj(Q)
  tng <- vapply(seq_len(ncol(cmp)), function(j) {
    dj <- d0 + eps * cmp[, j]
    dj <- dj - ray * sum(dj * ray); dj <- dj / sqrt(sum(dj^2))
    as.numeric(Pj(cbind(ray, dj)) - P0) / eps
  }, numeric(p * p))
  sv_t <- svd(tng)$d
  cat(sprintf("  tangent singular values: %s\n",
              paste(format(sv_t, digits = 3), collapse = "  ")))   # expect 1.41 1.41
  cat(sprintf("  -> measured PLANE dimension = %d (rotations = 3; the difference is\n",
              sum(sv_t > 1e-8 * max(sv_t))))
  cat("     the in-plane spin, 3 = 2 + 1). This is the correction to §5.1.\n")
  stopifnot(sum(sv_t > 1e-8 * max(sv_t)) == (p - 2L))
}

# ---- 5) What the arbitrary completion costs, WITHIN one pair ---------------
# The note's §5.4 measures the consequence ACROSS pairs. The tighter claim -- that
# the completion alone moves the answer at one fixed pair -- comes free: the
# package called twice on pair (1,2) gives two different distances.
cat("\n=== 5) Within-pair cost of the arbitrary completion (p = 4) ===\n")
{
  l12 <- bl_find_local_cf(S4$bl_res, S4$tgt, max_pairs = 1L, verbose = FALSE)
  lop <- bl_find_local_cf(S4$bl_res, S4$tgt, verbose = FALSE)
  a <- sqrt(l12$dist_mahalanobis)
  b <- sqrt(as.numeric(lop$all_distances_mahalanobis[["(1,2)"]]))
  cat(sprintf("  max_pairs = 1L         sqrt(mah) = %.6f\n", a))   # expect 1.096767
  cat(sprintf("  inside the full search sqrt(mah) = %.6f\n", b))   # expect 1.099395
  cat(sprintf("  difference %.2e -- same pair, same data, two calls.\n", abs(a - b)))
  cat("  REPORTED, not asserted: the size is data-dependent; the point is that it\n")
  cat("  is nonzero at all.\n")
}

# ---- 6) The in-plane objective (§6.3, stage 1) -----------------------------
# Nearest in-plane boundary point, by casting K rays from the target and bisecting
# each to its FIRST cutoff crossing. Correct for a nearest-point search: the
# nearest point is reached by some ray, and along that ray it must be the first
# crossing, since an earlier crossing would be nearer.
#
# VECTORISED ON PURPOSE. Every ray is evaluated at one step size in ONE predict
# call. A scalar version of this made 1,368 separate predict.gam calls per plane
# and was unusable -- a 600 s run had to be killed. Do not "simplify" it back.
ray_scan <- function(S, Vr, K = 36L, tmax = 6, march = 24L, nb = 14L) {
  tc  <- as.integer(S$pfun(matrix(S$x_t, 1L)) >= S$cutoff)
  z_t <- as.numeric(t(Vr) %*% S$x_t)
  a   <- seq(0, pi - pi / K, length.out = K)
  D   <- rbind(cbind(cos(a), sin(a)), -cbind(cos(a), sin(a)))    # 2K x 2
  cls_at <- function(s_vec)
    as.integer(S$pfun(sweep(D * s_vec, 2L, z_t, "+") %*% t(Vr)) >= S$cutoff)

  lo <- rep(0, nrow(D)); hi <- rep(NA_real_, nrow(D))
  for (s in seq(tmax / march, tmax, length.out = march)) {        # march to 1st crossing
    act <- is.na(hi); if (!any(act)) break
    cc <- cls_at(ifelse(act, s, lo)); nw <- act & (cc != tc)
    hi[nw] <- s; lo[act & !nw] <- s
  }
  ok <- !is.na(hi)
  if (!any(ok)) return(list(d2 = Inf, B_x = NULL))
  for (b in seq_len(nb)) {                                       # bisect
    mid <- (lo + hi) / 2
    cc  <- cls_at(ifelse(ok, mid, lo)); fl <- ok & (cc != tc)
    hi[fl] <- mid[fl]; lo[ok & !fl] <- mid[ok & !fl]
  }
  Xb <- sweep(D * hi, 2L, z_t, "+")[ok, , drop = FALSE] %*% t(Vr)
  d2 <- apply(Xb, 1L, function(xb)
          as.numeric(t(xb - S$x_t) %*% S$Minv %*% (xb - S$x_t)))  # SQUARED Mahalanobis
  k  <- which.min(d2)
  list(d2 = d2[k], B_x = as.numeric(Xb[k, ]))
}

# Deterministic orthonormal basis of the target's orthogonal complement, via a
# Householder reflector. No SVD, no sign ambiguity beyond the stated convention.
perp_basis <- function(x_t) {
  p  <- length(x_t)
  q  <- x_t / sqrt(sum(x_t^2))
  vv <- q - c(1, rep(0, p - 1))
  hh <- if (sum(vv^2) < 1e-24) diag(p) else diag(p) - 2 * outer(vv, vv) / sum(vv^2)
  U  <- hh[, -1, drop = FALSE]
  stopifnot(max(abs(t(U) %*% q)) < 1e-10,
            max(abs(crossprod(U) - diag(p - 1))) < 1e-10)
  list(q = q, U = U)
}

# span(q, U d), orthonormalised -- one admissible plane.
plane_from <- function(pb, d) {
  Vr <- cbind(pb$q, as.numeric(pb$U %*% d))
  Vr[, 2] <- Vr[, 2] - Vr[, 1] * sum(Vr[, 1] * Vr[, 2])
  Vr[, 2] <- Vr[, 2] / sqrt(sum(Vr[, 2]^2))
  Vr
}

# ---- 7) p = 3: the family is ONE parameter, so sweep it exhaustively -------
# Ranking settings (K = 18, march 12, bisect 8) -- coarser than section 8's,
# because stage 1 only has to rank planes; §6.3's stage 2 re-runs the package's
# own pipeline on the winner.
cat("\n=== 7) p = 3 exhaustive 1-degree sweep of the admissible pencil ===\n")
{
  pb3  <- perp_basis(S3$x_t)
  phis <- seq(0, pi - pi / 180, length.out = 180L)
  t0   <- proc.time()[["elapsed"]]
  dsw  <- vapply(phis, function(ph)
            ray_scan(S3, plane_from(pb3, c(cos(ph), sin(ph))),
                     K = 18L, march = 12L, nb = 8L)$d2, numeric(1))
  t_naive_sweep <- proc.time()[["elapsed"]] - t0

  prs3 <- gen_pairs(S3$p)
  lab  <- vapply(seq_len(nrow(prs3)), function(k)
            ray_scan(S3, min_rotation(S3$x_t, S3$V, prs3[k, ])$Vr,
                     K = 18L, march = 12L, nb = 8L)$d2, numeric(1))
  names(lab) <- apply(prs3, 1L, function(r) sprintf("(%d,%d)", r[1], r[2]))
  pk3 <- pkg_ref3                       # captured in section 0, before any pin

  cat(sprintf("  sweep    min %.6f | median %.6f | max %.6f\n",
              min(dsw), median(dsw), max(dsw)))   # expect 1.237877 / 1.490022 / 1.579964
  cat(sprintf("  argmin at phi = %.1f deg\n", phis[which.min(dsw)] * 180 / pi))  # expect 134.0
  cat(sprintf("  minimal-rotation labels: %s\n",
              paste(sprintf("%s %.6f", names(lab), lab), collapse = " | ")))
              # expect (1,2) 1.580762 | (1,3) 1.293857 | (2,3) 1.568419
  cat(sprintf("  the package, as it stands: %.6f\n", pk3))   # expect 1.312277

  # MONOTONICITY. Every plane the current method draws contains the target, so it
  # lies in the family -- the sweep's optimum cannot be worse than either.
  stopifnot(min(dsw) <= pk3 + 1e-9, min(dsw) <= min(lab) + 1e-9)
  cat(sprintf("  ASSERTED monotonicity: sweep <= package and <= best label.\n"))
  cat(sprintf("  improvement over the package: %.1f%% squared, %.1f%% distance\n",
              100 * (1 - min(dsw) / pk3), 100 * (1 - sqrt(min(dsw)) / sqrt(pk3))))

  nm <- sum(vapply(seq_along(dsw), function(k) {
    a <- dsw[if (k == 1L) 180L else k - 1L]
    b <- dsw[if (k == 180L) 1L else k + 1L]
    dsw[k] < a && dsw[k] < b }, logical(1)))
  cat(sprintf("  local minima on the curve: %d -> near-ties are real\n", nm))  # expect 2
  cat(sprintf("  pinning to (1,2) gives %.6f, within %.2f%% of the sweep's MAXIMUM\n",
              lab[["(1,2)"]], 100 * abs(lab[["(1,2)"]] / max(dsw) - 1)))
  cat("  -> a pinned pair is not a neutral default.\n")
  assign("t_naive_sweep", t_naive_sweep, envir = globalenv())
  assign("dsw", dsw, envir = globalenv())
}

# ---- 8) p = 4: two parameters, so a deterministic lattice ------------------
# Fibonacci/spiral hemisphere -- deterministic in N, no RNG anywhere.
cat("\n=== 8) p = 4 family lattice (200 points) ===\n")
{
  pb4  <- perp_basis(S4$x_t)
  lab4 <- vapply(seq_len(nrow(prs4)), function(k)
            ray_scan(S4, min_rotation(S4$x_t, S4$V, prs4[k, ])$Vr)$d2, numeric(1))
  names(lab4) <- pn4

  N  <- 200L; gr <- (1 + sqrt(5)) / 2
  zz <- (seq_len(N) - 0.5) / N                       # upper hemisphere only
  dd <- cbind(sqrt(1 - zz^2) * cos(2 * pi * seq_len(N) / gr),
              sqrt(1 - zz^2) * sin(2 * pi * seq_len(N) / gr), zz)
  fam <- apply(dd, 1L, function(d) ray_scan(S4, plane_from(pb4, d))$d2)
  pk4 <- pkg_ref4                       # captured in section 0, before any pin

  cat(sprintf("  6 minimal-rotation labels: %s\n",
              paste(sprintf("%s %.6f", names(lab4), lab4), collapse = " | ")))
  cat(sprintf("  best label %.6f | family %.6f | package %.6f\n",
              min(lab4), min(fam), pk4))
              # expect 1.147994 | 1.097519 | 1.147879
  stopifnot(min(fam) <= pk4 + 1e-9, min(fam) <= min(lab4) + 1e-9)
  cat(sprintf("  ASSERTED monotonicity. Improvement over best label: %.1f%%\n",
              100 * (1 - min(fam) / min(lab4))))
}

# ---- 8b) FINDING: the full search depends on WHAT WAS CALLED BEFORE IT ------
# Found while rebuilding this harness (2026-09-10) and not previously recorded.
# For one fixed target, bl_find_local_cf()'s reported distance changes according
# to the call history of the session -- specifically, a full search run
# immediately after a max_pairs = 1L call on the SAME objects returns a different
# distance from the same call made in isolation.
#
# This is NOT the floating-point under-determination of §4.1, and not §5.3's
# within-session pair flip: the selected pair is unchanged, 200 unrelated svd()
# calls have no effect, and building a second bl_result has no effect. It is
# specific to prior bl_find_local_cf() calls. The mechanism is not established
# here -- the installed package ships as a lazy-load database, so the source was
# not available to inspect -- so this is a measurement and a bug report, not a
# diagnosis.
cat("\n=== 8b) Order dependence of the full search (p = 4, row 49) ===\n")
{
  full <- function() bl_find_local_cf(S4$bl_res, S4$tgt, verbose = FALSE)
  pin  <- function() bl_find_local_cf(S4$bl_res, S4$tgt, max_pairs = 1L,
                                      verbose = FALSE)
  seqn <- list("full" = full, "pin" = pin, "full" = full, "full" = full,
               "pin" = pin, "full" = full, "full" = full)
  seen <- c()
  for (k in seq_along(seqn)) {
    lo <- seqn[[k]]()
    pr <- as.integer(lo$best_pair)
    cat(sprintf("  %d. %-5s -> pair (%d,%d) | dist_mah %.6f\n",
                k, names(seqn)[k], pr[1], pr[2], lo$dist_mahalanobis))
    if (names(seqn)[k] == "full") seen <- c(seen, lo$dist_mahalanobis)
  }
  # The SEQUENCE is itself not stable run to run -- two runs on 2026-09-10 gave
  # ... / 1.153789 / 1.147879 / ... and ... / 1.153789 / 1.153789 / ... at calls
  # 3-4. What is stable is the SET: the full search returns 1.147879 or 1.153789,
  # and the pinned search 1.208670 or 1.202898. Do not encode an expected order.
  cat(sprintf("  distinct full-search values in one session: %s\n",
              paste(unique(sprintf("%.6f", seen)), collapse = " | ")))
  if (length(unique(sprintf("%.6f", seen))) > 1L)
    cat("  -> the SAME call on the SAME objects returns more than one answer.\n")
  cat("  REPORTED, not asserted: this is the behaviour under test, so asserting a\n")
  cat("  particular sequence would freeze in a bug rather than surface it.\n")
  cat("  Practical rule: capture a quoted counterfactual from a FIRST call on a\n")
  cat("  fresh result object, and re-run to confirm it is stable (§5.3, §6.5 test 2).\n")
}

# ---- 9) Cost -- REPORTED ONLY ---------------------------------------------
# §6.3's cost paragraph. The plane-at-a-time sweep is dominated by per-call
# overhead, not arithmetic: the package's grid does 40,000 rows in one vectorised
# call per pair, while a per-plane ray scan does 36 rows per call. Batching the
# scan ACROSS planes puts the whole sweep in one call per march/bisect step.
cat("\n=== 9) Cost (reported, not asserted -- wall clock is machine-dependent) ===\n")
{
  t0 <- proc.time()[["elapsed"]]
  invisible(bl_find_local_cf(S3$bl_res, S3$tgt, verbose = FALSE))
  t_pkg <- proc.time()[["elapsed"]] - t0

  # Batched: every plane's every ray in ONE predict call per step.
  pb3 <- perp_basis(S3$x_t)
  NP <- 180L; K <- 18L; march <- 12L; nb <- 8L; tmax <- 6
  phis <- seq(0, pi - pi / NP, length.out = NP)
  VRs <- lapply(phis, function(ph) plane_from(pb3, c(cos(ph), sin(ph))))
  aa  <- seq(0, pi - pi / K, length.out = K)
  D   <- rbind(cbind(cos(aa), sin(aa)), -cbind(cos(aa), sin(aa)))
  nR  <- nrow(D)
  Zt  <- t(vapply(VRs, function(Vr) as.numeric(t(Vr) %*% S3$x_t), numeric(2)))
  pl_idx <- rep(seq_len(NP), each = nR); ray_idx <- rep(seq_len(nR), times = NP)
  tc <- as.integer(S3$pfun(matrix(S3$x_t, 1L)) >= S3$cutoff)

  cls_all <- function(s_vec) {
    Zp <- Zt[pl_idx, , drop = FALSE] + D[ray_idx, , drop = FALSE] * s_vec
    Xp <- matrix(0, length(s_vec), S3$p)
    for (k in seq_len(NP)) {
      m <- pl_idx == k
      Xp[m, ] <- Zp[m, , drop = FALSE] %*% t(VRs[[k]])
    }
    as.integer(S3$pfun(Xp) >= S3$cutoff)
  }

  t0 <- proc.time()[["elapsed"]]
  lo <- rep(0, NP * nR); hi <- rep(NA_real_, NP * nR)
  for (s in seq(tmax / march, tmax, length.out = march)) {
    act <- is.na(hi); if (!any(act)) break
    cc <- cls_all(ifelse(act, s, lo)); nw <- act & (cc != tc)
    hi[nw] <- s; lo[act & !nw] <- s
  }
  ok <- !is.na(hi)
  for (b in seq_len(nb)) {
    mid <- (lo + hi) / 2
    cc <- cls_all(ifelse(ok, mid, lo)); fl <- ok & (cc != tc)
    hi[fl] <- mid[fl]; lo[ok & !fl] <- mid[ok & !fl]
  }
  Zf <- Zt[pl_idx, , drop = FALSE] + D[ray_idx, , drop = FALSE] * hi
  best <- rep(Inf, NP)
  for (k in seq_len(NP)) {
    m <- (pl_idx == k) & ok
    if (!any(m)) next
    Xb <- Zf[m, , drop = FALSE] %*% t(VRs[[k]])
    best[k] <- min(apply(Xb, 1L, function(xb)
      as.numeric(t(xb - S3$x_t) %*% S3$Minv %*% (xb - S3$x_t))))
  }
  t_batched <- proc.time()[["elapsed"]] - t0

  cat(sprintf("  package full search (3 pairs, 200x200 grid) : %6.2f s\n", t_pkg))
  cat(sprintf("  sweep, plane at a time                     : %6.2f s  (%.1fx)\n",
              t_naive_sweep, t_naive_sweep / t_pkg))          # recorded 6.3x
  cat(sprintf("  sweep, batched across planes               : %6.2f s  (%.2fx)\n",
              t_batched, t_batched / t_pkg))                  # recorded 1.65x
  cat(sprintf("  batched rows per call: %d planes x %d rays = %d, in %d calls\n",
              NP, nR, NP * nR, march + nb))
  cat(sprintf("  ray-scan predictions per plane at K = 36: %d\n", 2L * 36L * (24L + 14L)))
  cat(sprintf("  batched sweep agrees with section 7: min %.6f\n", min(best)))
  stopifnot(abs(min(best) - min(dsw)) < 1e-9)   # same objective, different batching
  cat(sprintf(paste0("  -> the search is a real COST, about %.1fx batched",
                     " (%.0fx unbatched). Not a saving.\n"),
              t_batched / t_pkg, t_naive_sweep / t_pkg))
  cat("     Recorded on 2026-09-04: 1.65x batched, 6.3x unbatched. Both figures\n")
  cat("     move with machine load -- the ratio, not the seconds, is the point.\n")
}

# ---- 10) The thesis's worked example is preserved (§6.2) -------------------
# At p = 2, r = 1 the admissible family is a single line -- there is exactly one
# line through the origin containing the target -- so the minimal rotation and
# Procrustes MUST agree, and the §2.5 exposition is untouched by the change.
#
# The five points are 17_1a's (originally A3's), defined inline rather than
# sourced: 17_1a draws its six-panel figure on load.
cat("\n=== 10) Thesis reproduction (p = 2, r = 1, target row 5) ===\n")
{
  X2   <- matrix(c(1.5, 1.5,  4, 2,  2, 4,  4, 4,  3, 4), ncol = 2, byrow = TRUE)
  X2c  <- scale(X2, center = TRUE, scale = FALSE)     # centred, NOT scaled
  V2   <- svd(X2c)$v
  xt2  <- as.numeric(X2c[5L, ])

  Y  <- rbind(-xt2, rep(0, 2), xt2)
  Z  <- cbind(Y %*% V2[, 1, drop = FALSE], matrix(0, 3L, 1L))
  sv <- svd(t(Y %*% V2) %*% Z)
  th_proc <- acos(pmin(1, abs((sv$v %*% t(sv$u))[1, 1]))) * 180 / pi

  xP <- V2[, 1] * sum(xt2 * V2[, 1])
  th_min <- acos(min(1, sqrt(sum(xP^2)) / sqrt(sum(xt2^2)))) * 180 / pi

  cat(sprintf("  Procrustes       %.4f deg\n", th_proc))   # expect 28.1411
  cat(sprintf("  minimal rotation %.4f deg\n", th_min))    # expect 28.1411
  cat(sprintf("  thesis (§2.5)    28.14 deg | agreement %.2e\n",
              abs(th_min - th_proc)))                      # expect 6.75e-14
  stopifnot(abs(th_min - th_proc) < 1e-8, abs(th_min - 28.14) < 0.01)
  cat("  ASSERTED: the two constructions agree, so the change removes only the\n")
  cat("  freedom that exists from p = 3 upward.\n")
}

cat("\n=== Done. Asserted: sections 2, 3, 4, 7, 8, 10. Reported: 5, 9. ===\n")
cat("Still unrun (see §6.5): test 5 (variance retained), test 7 at p >= 4,\n")
cat("and cross-BLAS determinism.\n")
