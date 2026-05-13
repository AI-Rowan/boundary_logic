############################################################
# Dataset documentation
# Roxygen2 entries for files bundled in inst/extdata/
############################################################

#' Loan Default Dataset
#'
#' A dataset of 45,000 personal loan applications used to demonstrate
#' the \pkg{boundarylogic} workflow on a credit-risk classification problem.
#'
#' @format A data frame with 45,000 rows and 14 variables:
#' \describe{
#'   \item{person_age}{Age of the applicant (numeric).}
#'   \item{person_gender}{Gender of the applicant: \code{"female"} or
#'     \code{"male"}. Integer-encoded in the example script as 0 = female,
#'     1 = male.}
#'   \item{person_education}{Highest education level attained: \code{"High School"},
#'     \code{"Associate"}, \code{"Bachelor"}, \code{"Master"}, or
#'     \code{"Doctorate"}. Integer-encoded as an ordinal 0–4.}
#'   \item{person_income}{Annual income in USD (numeric).}
#'   \item{person_emp_exp}{Years of employment experience (numeric).}
#'   \item{person_home_ownership}{Home ownership status: \code{"MORTGAGE"},
#'     \code{"OTHER"}, \code{"OWN"}, or \code{"RENT"}. Integer-encoded as
#'     0–3 in the example script.}
#'   \item{loan_amnt}{Loan amount applied for in USD (numeric).}
#'   \item{loan_intent}{Stated purpose of the loan: \code{"DEBTCONSOLIDATION"},
#'     \code{"EDUCATION"}, \code{"HOMEIMPROVEMENT"}, \code{"MEDICAL"},
#'     \code{"PERSONAL"}, or \code{"VENTURE"}. Integer-encoded as 0–5.}
#'   \item{loan_int_rate}{Interest rate of the loan as a percentage (numeric).}
#'   \item{loan_percent_income}{Loan amount as a fraction of annual income
#'     (numeric).}
#'   \item{cb_person_cred_hist_length}{Length of credit history in years
#'     (numeric).}
#'   \item{credit_score}{Applicant credit score (numeric).}
#'   \item{previous_loan_defaults_on_file}{Whether a prior default is on record:
#'     \code{"No"} or \code{"Yes"}. Integer-encoded as 0 = No, 1 = Yes.}
#'   \item{loan_status}{Binary outcome: 1 = default, 0 = no default.}
#' }
#'
#' @details Access the raw file via:
#' \preformatted{
#' read.csv(system.file("extdata", "loan_data.csv", package = "boundarylogic"))
#' }
#'
#' The raw file contains character columns for \code{person_gender},
#' \code{person_education}, \code{person_home_ownership}, \code{loan_intent},
#' and \code{previous_loan_defaults_on_file}. These must be integer-encoded
#' before passing to \code{bl_prepare_data()}. The example script
#' \code{scripts/03_loan_status_Boundary_Logic.R} shows the encoding step
#' and an iterative three-phase workflow, including a domain filter that
#' retains only applicants with a prior default on file.
#'
#' @source \url{https://github.com/TSMathi/loan_approval_analysis/tree/main}
#'
#' @name loan_data
#' @docType data
#' @keywords datasets
NULL
