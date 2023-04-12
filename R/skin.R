#' @title Skin Irritation Data
#'
#' @description The data set is a re-simulated part from an ongoing
#' neuodermatitis study where researchers investigate the
#' efficacy of an ointment in reducing the severity of skin irritation on the
#' backs of the hands of 25 neurodermatitis patients, where 10 patients’ backs
#' of the hands were rubbed with the ointment and 15 were not. The response is
#' a BI-RADS rating score and the lower the score the better the clinical record.
#' Every remarkable skin irritation was graded on every patients back of the
#' hands and thus, the numbers of replicates differ across the patients.
#'
#' @usage skin
#'
#' @format A data frame with 107 rows and 4 variables:
#' \describe{
#'   \item{\code{tx }}{treatment group}
#'   \item{\code{intervention }}{intervention period indicator}
#'   \item{\code{subject }}{subject ID}
#'   \item{\code{score }}{BI-RADS rating score, where 1 = very mild irritation,
#'   2 =slight irritation, 3 =mild irritation, 4 =heavy irritation and
#' 5 =severe irritation}
#' }
#'
#' @references Roy, A, Harrar, SW, Konietschke, F. The nonparametric
#' Behrens-Fisher problem with dependent replicates. Statistics in Medicine.
#' 2019; 38: 4939– 4962. https://doi.org/10.1002/sim.8343
#'
#' @examples
#' data(skin)
#' head(skin)
"skin"

