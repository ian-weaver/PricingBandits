#' @keywords internal
"_PACKAGE"

#' @import dplyr
#' @import Matrix
#' @importFrom stats na.omit rbeta runif sd
#' @importFrom nloptr nloptr
#' @importFrom MASS mvrnorm
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom R.utils withTimeout
#' @importFrom hash hash
NULL

# quiet R CMD check NOTEs about dplyr non-standard evaluation columns
utils::globalVariables(c("s_t", "n_t", "s_t.y", "n_t.y", "PricesTested", "PurchaseDecisions"))
