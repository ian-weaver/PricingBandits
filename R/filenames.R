#' @title Extract the Name of a Policy Function
#' @description This function takes a policy function as input and returns its name as a string.
#' @param Policy A function representing the policy whose name is to be extracted.
#' @return A character vector containing the name of the policy function.
#' @examples
#' MyPolicy <- UCB
#' PolicyAsString(MyPolicy)
#' @export
PolicyAsString <- function(Policy){
  env_vars <- ls(environment(Policy))
  env_vars <- setdiff(env_vars, "Policy")  # Remove 'Policy' itself
  function_names <- env_vars[sapply(env_vars, function(x) identical(get(x, environment(Policy)), Policy))]
  return (function_names)
}
