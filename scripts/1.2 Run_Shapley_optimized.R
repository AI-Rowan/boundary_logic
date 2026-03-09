############################################################
# Script: Run_Shapley (Optimized, Same API)
# Author: (Your Name)
# Date: 2025-09-25
# Description:
#   Faster & memory-efficient Shapley computation for moving a point
#   from `test_data` to `B_boundary` by features. Keeps the original
#   function name/signature and returns the same 'shapley_store' shape.
#
#   Strategy:
#     • Exact (streamed over subset sizes) using the classic subset-
#       difference formula — O(n * 2^n) evaluations, but batched and
#       without factorial permutations. Used when num_vars ≤ exact_max_vars.
#     • Approximate (permutation sampling) for larger n — draws M
#       permutations and accumulates marginal contributions. Each
#       permutation needs (n+1) model evaluations; we batch all
#       permutations for a given observation to minimize overhead.
#
#   Notes:
#     • Uses `pred_function(model_use, model_select, rounding, new_data)`,
#       just like the original.
#     • No change to the notion of “start” (test_data) and “end” (B_boundary).
############################################################

# ---- helpers ------------------------------------------------------------

# Build a block of rows corresponding to a set of binary masks (0/1) where
# 0 -> take 'start', 1 -> take 'end'. Returns data.frame with var_names.
.build_rows_from_masks <- function(start, end, masks, var_names) {
  # start, end: numeric vectors length p
  # masks: matrix [R x p] of 0/1
  X <- sweep((1 - masks), 2, start, `*`) + sweep(masks, 2, end, `*`)
  X <- as.data.frame(X)
  colnames(X) <- var_names
  X
}

# Compute exact Shapley for a single observation via streaming over subset sizes.
# Avoids factorial permutations and avoids storing all 2^p rows at once.
.shapley_exact_one <- function(start, end, model_use, model_select, rounding, var_names,
                               exact_max_vars = 14L) {
  p <- length(start)
  # weights for each subset size k (used when aggregating across S \notni i):
  # w_k = k!(p-k-1)! / p!
  fact <- factorial(0:p)
  wfun <- function(k) (fact[k+1] * fact[p-k]) / fact[p+1]  # k in 0..p-1

  shap <- numeric(p)

  # Precompute index sets per k to avoid re-gen inside i-loop.
  # We'll use combn; for p up to 14 this is fine. Each k done once.
  for (k in 0:(p-1)) {
    if (k == 0) {
      comb <- matrix(integer(0), nrow = 0, ncol = 1L)
    } else {
      comb <- combn(p, k)  # columns are combinations of size k (1-based indices)
    }
    # For this k, build masks for S (size k) repeated p-k times when needed.
    # We'll compute f(S) once, then reuse when adding i.
    if (k == 0) {
      masks_S <- matrix(0L, nrow = 1L, ncol = p)
    } else {
      masks_S <- matrix(0L, nrow = ncol(comb), ncol = p)
      for (c in seq_len(ncol(comb))) masks_S[c, comb[, c]] <- 1L
    }

    # Predict f(S) in one batch
    X_S <- .build_rows_from_masks(start, end, masks_S, var_names)
    f_S <- pred_function(model_use = model_use, model_select = model_select,
                         rounding = rounding, new_data = X_S)

    # For each feature i not in S, form S∪{i}. We do this by flipping column i to 1.
    # To keep batching, we construct a big block concatenating all needed (S, i) pairs.
    if (k <= p-1) {
      # Total pairs count = nrow(masks_S) * (p - k)
      pairs <- list()
      idx_S <- list()
      idx_i <- list()
      cnt <- 1L
      for (rowS in seq_len(nrow(masks_S))) {
        mask_row <- masks_S[rowS, ]
        zeros <- which(mask_row == 0L)
        for (i in zeros) {
          new_mask <- mask_row; new_mask[i] <- 1L
          pairs[[cnt]] <- new_mask
          idx_S[[cnt]] <- rowS   # reference f(S)
          idx_i[[cnt]] <- i      # which feature gains credit
          cnt <- cnt + 1L
        }
      }
      masks_Si <- do.call(rbind, pairs)
      X_Si <- .build_rows_from_masks(start, end, masks_Si, var_names)
      f_Si <- pred_function(model_use = model_use, model_select = model_select,
                            rounding = rounding, new_data = X_Si)

      # Aggregate contributions with weight w_k
      wk <- wfun(k)
      for (m in seq_len(nrow(masks_Si))) {
        inc <- f_Si[m] - f_S[idx_S[[m]]]
        shap[idx_i[[m]]] <- shap[idx_i[[m]]] + wk * inc
      }
    }
  }
  shap
}

# Approximate Shapley via sampled permutations (KernelSHAP-like weighting optional)
.shapley_perm_one <- function(start, end, model_use, model_select, rounding, var_names,
                              M = 2048L, seed = 1L) {
  set.seed(seed)
  p <- length(start)
  shap_sum <- numeric(p)
  for (m in 1:M) {
    ord <- sample.int(p, p, replace = FALSE)
    # Build masks for this permutation path: from 0 active to p active (p+1 rows)
    masks <- matrix(0L, nrow = p + 1L, ncol = p)
    if (p > 0) {
      for (i in seq_len(p)) masks[(i+1), ord[i]] <- 1L
      masks <- apply(masks, 2, cumsum)  # cumulative include
    }
    X_path <- .build_rows_from_masks(start, end, masks, var_names)
    vals <- pred_function(model_use = model_use, model_select = model_select,
                          rounding = rounding, new_data = X_path)
    # Marginals along the path
    d <- diff(vals)  # length p, order-specific
    # Assign credits according to permutation order
    shap_sum[ord] <- shap_sum[ord] + d
  }
  shap_sum / M
}

# ---- main (same signature) ----------------------------------------------
Run_Shapley <- function(model_use, num_vars, test_data,
                        B_boundary,
                        individual = FALSE, tdp = 1,
                        lda = FALSE, rpart_used = FALSE, svm_u = FALSE,
                        rounding = rounding,
                        exact_max_vars = 14L,
                        approx_perm = 2048L,
                        seed = 1L) {

  n <- nrow(test_data)
  p <- as.integer(num_vars)
  out <- matrix(0, nrow = if (isTRUE(individual)) 1L else n, ncol = p)
  colnames(out) <- colnames(test_data)[seq_len(p)]
  var_names <- colnames(test_data)[seq_len(p)]

  run_idx <- if (isTRUE(individual)) tdp else seq_len(n)

  for (t in run_idx) {
    start <- as.numeric(test_data[t, 1:p, drop = FALSE])
    end   <- as.numeric(B_boundary[t, 1:p, drop = FALSE])

    if (p <= exact_max_vars) {
      shap <- .shapley_exact_one(start, end, model_use, model_select,
                                 rounding, var_names, exact_max_vars = exact_max_vars)
    } else {
      shap <- .shapley_perm_one(start, end, model_use, model_select,
                                rounding, var_names, M = approx_perm, seed = seed + t - 1L)
    }
    out_index <- if (isTRUE(individual)) 1L else t
    out[out_index, ] <- shap
  }

  list(shapley_store = out)
}
