# ============================================================================
# bl_select_target(target = i, data = df): which observation's class colours
# the target?
# ----------------------------------------------------------------------------
# Reproduces every number and the figure in select_target_class_note.md.
#
# Set-up: the Paper 1 Pima set-up -- 392 complete rows, p = 7 (BloodPressure
# dropped), 80/20 split with seed 121, no outlier trimming, GAM in raw units,
# CVA on the four prediction outcomes.
#
# Run from this folder:
#   Rscript select_target_class_example.R
# To measure the working tree rather than the installed build, run
# devtools::load_all("<package root>") first and source() this file with
# chdir = TRUE.
#
# Outputs (written next to this script):
#   select_target_class_results.txt   the counts and the chosen target
#   select_target_class_figure.png    the two-panel figure
# ============================================================================

if (!"package:boundarylogic" %in% search()) library(boundarylogic)

# ---- 1. Data, model, result (Paper 1 set-up) --------------------------------
pima_raw <- read.csv(system.file("extdata", "pima_diabetes.csv",
                                 package = "boundarylogic"))
keep <- with(pima_raw, Insulin != 0 & SkinThickness != 0 &
                       Glucose != 0 & BMI != 0)
pima <- pima_raw[keep, ]
stopifnot(nrow(pima) == 392L)

feature_cols <- c("Pregnancies", "Glucose", "SkinThickness", "Insulin",
                  "BMI", "DiabetesPedigreeFunction", "Age")

bl_dat <- bl_prepare_data(pima, class_col = "Outcome",
                          feature_cols = feature_cols,
                          train_fraction = 0.8, seed = 121L,
                          hull_fraction = 1, verbose = FALSE)

gam_fit <- mgcv::gam(
  class ~ s(Pregnancies) + s(Glucose) + s(SkinThickness) + s(Insulin) +
    s(BMI) + s(DiabetesPedigreeFunction) + s(Age),
  family = binomial, data = bl_dat$train_data)

bl_mod <- bl_wrap_model(gam_fit, model_type = "GAM",
                        var_names = bl_dat$var_names,
                        train_data = bl_dat$train_data)

bl_results <- bl_build_result(bl_dat, bl_mod, method = "CVA")

train_data <- bl_results$train_data
test_data  <- bl_results$test_data
n_train    <- nrow(train_data)
n_test     <- nrow(test_data)

# ---- 2. Class mismatch: train_data[i] against test_data[i] -----------------
# The local-biplot colour is the confusion group of (true class, predicted
# class). The predicted class is the target's own, so the colour is wrong
# exactly when the class read from test_data[i] differs from train_data[i]'s.
# For i > n_test the test_data lookup is out of range, so the target loses its
# confusion colour and is drawn by its predicted class only.
idx_in    <- seq_len(min(n_train, n_test))
mismatch  <- train_data$class[idx_in] != test_data$class[idx_in]
n_mis     <- sum(mismatch)
n_out     <- max(n_train - n_test, 0L)

# ---- 3. One target from train_data with a mismatched class -----------------
confusion_col <- function(true_class, pred_class) {
  if (is.na(true_class)) return(if (pred_class == 1L) "red" else "blue")
  if (true_class == 1L && pred_class == 1L) return("red")
  if (true_class == 0L && pred_class == 0L) return("blue")
  if (true_class == 0L && pred_class == 1L) return("purple")
  "orange"
}
confusion_lab <- c(red = "TP", blue = "TN", purple = "FP", orange = "FN")

# First mismatched row whose local search finds a counterfactual.
target_i <- NA_integer_
for (i in which(mismatch)) {
  tgt      <- bl_select_target(bl_results, target = i, data = train_data)
  bl_local <- bl_find_local_cf(bl_results, tgt, verbose = FALSE)
  if (isTRUE(bl_local$solution_found)) { target_i <- i; break }
}
stopifnot(!is.na(target_i))

class_own  <- as.integer(train_data$class[target_i])
class_test <- as.integer(test_data$class[target_i])
pred_class <- tgt$pred_class
col_wrong  <- confusion_col(class_test, pred_class)
col_right  <- confusion_col(class_own,  pred_class)

# ---- 4. Figure: the colour from test_data[i] beside the target's own -------
# plot.bl_local_result() colours the target from bl_target$class_obs. Setting
# class_obs on a copy reproduces both colourings from one search:
#   left  -- test_data[i]'s class, the colour drawn when the class is looked up
#            in test_data by row number;
#   right -- train_data[i]'s class, the target's own.
local_wrong <- bl_local
local_wrong$bl_target$class_obs <- class_test
local_right <- bl_local
local_right$bl_target$class_obs <- class_own

png("select_target_class_figure.png", width = 1600L, height = 820L, res = 130L)
graphics::par(mfrow = c(1L, 2L))
plot(local_wrong, print_summary = FALSE,
     new_title = sprintf("Class from test_data[%d]: %s (%s)", target_i,
                         confusion_lab[[col_wrong]], col_wrong))
plot(local_right, print_summary = FALSE,
     new_title = sprintf("Class from train_data[%d]: %s (%s)", target_i,
                         confusion_lab[[col_right]], col_right))
grDevices::dev.off()

# ---- 5. Results -------------------------------------------------------------
results <- c(
  sprintf("boundarylogic version : %s", utils::packageVersion("boundarylogic")),
  sprintf("rows after zero filter: %d", nrow(pima)),
  sprintf("training rows         : %d", n_train),
  sprintf("test rows             : %d", n_test),
  sprintf("rows compared (i <= n_test)          : %d", length(idx_in)),
  sprintf("class differs, train[i] vs test[i]   : %d of %d", n_mis, length(idx_in)),
  sprintf("rows with i > n_test (lookup out of range, colour by prediction): %d",
          n_out),
  sprintf("target row (train_data)              : %d", target_i),
  sprintf("target pred_prob / pred_class        : %.3f / %d",
          tgt$pred_prob, pred_class),
  sprintf("class of train_data[i] / test_data[i]: %d / %d", class_own, class_test),
  sprintf("colour from test_data[i]             : %s (%s)",
          confusion_lab[[col_wrong]], col_wrong),
  sprintf("colour from the target's own class   : %s (%s)",
          confusion_lab[[col_right]], col_right),
  sprintf("bl_target has class_obs field        : %s",
          "class_obs" %in% names(tgt))
)
writeLines(results, "select_target_class_results.txt")
cat(results, sep = "\n")
