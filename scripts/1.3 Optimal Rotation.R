############################################################
# Section 1.3: Optimal Rotation for Biplot (V → Vrho)
# Description:
#   Computes a rotation matrix A that aligns 2D biplot plane 
#   spanned by `proj` (typically PCs 1 & 2) with a selected target point (tdp)
#   Returns the rotated loading matrices Vrho and tVrho for downstream use.
#
# Notes:

############################################################

###################################
####### Apply rho to projection 3d
###################################

#' Rotate the loading matrix so a target point lies on the biplot plane
#'
#' Constructs a 3-row matrix Y consisting of (-x_tdp, 0, +x_tdp) in standardized
#' feature space, then solves for a rotation A via SVD to map from full space
#' to the reduced (proj) space. Returns Vrho = V %*% t(A) and its inverse.
#'
#' @param X Data frame/matrix with features in the first num_vars columns
#' @param tdp Integer index of target data point for alignment
#' @param standardise Logical; whether to standardize using Xcenter/Xsd
#' @param Xcenter Numeric vector of feature means
#' @param Xsd Numeric vector of feature sds
#' @param CVA_method Logical; if TRUE forces standardise = FALSE
#' @param proj Integer vector of projection indices (e.g., c(1,2))
#' @param V Loading/eigenvector matrix (columns correspond to components)
#'
#' @return list(Vrho = rotated loadings, tVrho = inverse, A = rotation matrix)
Biplot_rotation <- function(X = biplot_data,
                            tdp = 1,
                            standardise = standardise_use, # set to false for CVA plots
                            Xcenter = X_center,
                            Xsd = X_sd,
                            CVA_method = CVA_method,
                            proj = pc.values,
                            V = V) {

  # Select scaling based on PCA/CVA choice
  sv <- Xsd
  if(CVA_method == T) standardise  <- F # set to false if CVA is true
  if(standardise  == F) sv = F

  # Standardize input (applies scaling to X; outputs scaled coordinates)
  X.st <- scale(X[, 1:num_vars], center = Xcenter, scale = sv)

  # Construct Y with three rows: -x_tdp, 0, +x_tdp
  Y <- as.matrix(rbind(-X.st[tdp,],
                       rep(0, num_vars),
                       X.st[tdp,]))

  # Extract loadings for full and reduced (proj) spaces
  V  <- V           # full dimensions
  Vr <- V[, proj]   # reduced biplot dimensions

  # Map Y to full PCA/CVA space and to reduced space
  YV  <- Y %*% V
  r   <- 2                         # target biplot dimensionality
  YVr <- Y %*% Vr                  # reduced-space coords
  if(num_vars > r) YVr <- cbind(YVr, matrix(0, ncol = num_vars - r, nrow = 3))

  # Solve for rotation A using SVD of cross-covariance
  svdyx <- svd(t(YV) %*% YVr)
  A     <- svdyx$v %*% t(svdyx$u)

  # Rho diagnostic (kept for reference with original code)
  rho <- sum(diag(YVr %*% A %*% t(YV))) / sum(diag(t(YVr) %*% (YVr)))

  # Apply rotation to V and compute inverse
  Vrho  <- V %*% t(A)           # rotated loading matrix
  tVrho <- solve(Vrho)

  # Optionally bypass PCA/CVA (identity) if requested upstream
  if(no_pca  == T) Vrho <- diag(num_vars)

  return(list(Vrho = Vrho,
              tVrho = tVrho,
              A = A))
}

#########################
####### END #############
#########################
