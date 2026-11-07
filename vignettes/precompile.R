# Precompute the loan-default vignette.
#
# The vignette runs the whole workflow on the full loan data, so it is knitted
# here, by the maintainer, and shipped with its output and figures:
#
#   vignettes/Boundary_Logic_loan_default_workflow.Rmd.orig  the source (edit this)
#   vignettes/Boundary_Logic_loan_default_workflow.Rmd       knitted, shipped
#   vignettes/figures/loan-*.png                             knitted figures
#
# R CMD build and pkgdown render the knitted .Rmd, which runs no code. Re-run
# this script after any change to the package code or to the .Rmd.orig, from
# the package root:
#
#   Rscript vignettes/precompile.R
#
# It loads the development copy of the package (never an installed one) and
# needs the optional packages fastshap and shapviz for the SHAP comparison.

local({
  if (!file.exists("DESCRIPTION") ||
      read.dcf("DESCRIPTION", fields = "Package")[1, 1] != "boundarylogic")
    stop("Run this script from the boundarylogic package root.", call. = FALSE)
  root <- normalizePath(".", winslash = "/")

  for (pkg in c("pkgload", "knitr", "fastshap", "shapviz"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is needed to precompute the vignette.",
           call. = FALSE)

  pkgload::load_all(root, quiet = TRUE)
  if (!pkgload::is_dev_package("boundarylogic"))
    stop("boundarylogic is not loaded from the development copy.", call. = FALSE)

  old_wd <- setwd(file.path(root, "vignettes"))
  on.exit(setwd(old_wd), add = TRUE)

  src <- "Boundary_Logic_loan_default_workflow.Rmd.orig"
  out <- "Boundary_Logic_loan_default_workflow.Rmd"

  # Remove the figures of the previous run, so none goes stale.
  unlink(Sys.glob("figures/loan-*.png"))

  t0 <- Sys.time()
  knitr::knit(src, output = out, envir = new.env(parent = globalenv()),
              quiet = TRUE)
  message(sprintf("Knitted %s in %.1f s.", out,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  # Checks: no live chunk is left, and every figure referenced exists.
  txt <- readLines(out, encoding = "UTF-8")
  if (any(grepl("^```[{]r", txt)))
    stop("The knitted vignette still has live R chunks.", call. = FALSE)
  figs <- unlist(regmatches(txt, gregexpr("figures/loan-[^)\"' ]+[.]png", txt)))
  missing <- figs[!file.exists(figs)]
  if (length(missing))
    stop("Missing figures: ", paste(missing, collapse = ", "), call. = FALSE)
  unused <- setdiff(Sys.glob("figures/loan-*.png"), figs)
  if (length(unused))
    warning("Figures not referenced: ", paste(unused, collapse = ", "),
            call. = FALSE)
  message(sprintf("%d figures referenced, all present.", length(unique(figs))))
})
