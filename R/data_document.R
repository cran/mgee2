#' obs1: simulated observed data
#'
#' @format a dataframe with 3000 rows and 8 variables
#'  \describe{
#'     \item{ID}{individual id number}
#'     \item{Y}{true response, factor variable}
#'     \item{X}{true error-prone covariate, factor variable}
#'     \item{treatment}{error-free covariate}
#'     \item{visit}{serial number of each visit}
#'     \item{S}{observed response, same as Y when in the validation set(delta=1)}
#'     \item{W}{observed error-prone covariate, same as X when in the validation
#'              set (delta=1)}
#'     \item{delta}{indicator variable, 1 if in the validation set, 0 if not.}
#'    }
"obs1"


