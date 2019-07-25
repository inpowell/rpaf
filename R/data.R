#' Mini-Finland Health Survey
#'
#' A dataset containing various health indicators and outcomes from the
#' Mini-Finland Health Survey cohort study by Aromaa et al. (1989).
#'
#' @format A data frame with 4517 rows and 14 variables:
#' \describe{
#' \item{ID}{unique identifier of each person in the study}
#' \item{DEATH}{\code{TRUE} if the subject died during follow-up}
#' \item{DEATH_FT}{length of time until death or censoring}
#' \item{DIAB}{\code{TRUE} if the subject contraced Type II Diabetes during follow-up}
#' \item{DIAB_FT}{length of time until diabetes occurrence or censoring}
#' \item{BYEAR}{year of birth}
#' \item{B_COHORT}{birth cohort, in blocks of 20 years from 1890 to 1930}
#' \item{AGE}{age in years at baseline}
#' \item{AGEGRP}{age group (years at baseline), in blocks of ten years from 50-79}
#' \item{SEX}{sex of participant (male/female)}
#' \item{BMI}{body mass index, in kg/m^2}
#' \item{BMI_2}{indicator of high BMI (>=25.0 kg/m^2) or normal BMI (<25.0 kg/m^2)}
#' \item{BP}{indicator for normal or elevated blood pressure}
#' \item{SMOKE}{smoking status: never smoked, former smoker, pipe/cigar only or
#'              <30 cigarettes/day, or >=30 cigarettes/day}
#' }
#'
#' @section References:
#'
#' \itemize{
#'
#' \item Aromaa A, Heliövaara M, Impivaara O, Knekt P, Maatela J (1989).
#'   “Aims, Methods and Study Population. Part 1.” In A Aromaa, M Heliövaara, O
#'   Impivaara, P Knekt, J Maatela (eds.), The Execution of the Mini-Finland
#'   Health Survey. (In Finnish, English summary). Publications of the Social
#'   Insurance Institution, Finland, Helsinki and Turku, ML:88.
#'
#' }
"minifhs"
