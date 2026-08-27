#' @title Basis Function
#' @description Computes the values of a basis function at a test point, given a number of knots.
#' @param x Test point where the basis function is evaluated.
#' @param J Number of knots for the basis function.
#' @return A vector of length J containing the basis function values at the test point.
#' @details The function computes the basis function values using a piecewise linear approach with J knots. The calculation adjusts based on the location of the test point relative to the knot positions.
#' @export
BasisFunction <- function(x, J){
  
  Delta_J = 1/(J-1)
  Knots   = c(0, seq(J-1)/(J-1))
  Output  = rep(0, J)
  i       = max(which(Knots <= x))
  
  if(i == 1){
    Output[1] = x - 0.5*(x^2)/Delta_J
    Output[2] = x - Knots[2]*x/Delta_J + 0.5*x^2/Delta_J
  }
  if(i == 2){
    Output[1] = Delta_J/2
    Output[2] = Delta_J/2 + (x-Knots[2])*(1+Knots[2]/Delta_J) - 0.5*(x^2-Knots[2]^2)/Delta_J
    Output[3] = (x-Knots[2])*(1-Knots[3]/Delta_J) + 0.5*(x^2-Knots[2]^2)/Delta_J
  }
  if(i == J){
    Output[1]       = Delta_J/2
    Output[2:(J-1)] = Delta_J
    Output[J]       = Delta_J/2
  }
  if(i != 1 && i != 2 && i != J){
    Output[1]       = Delta_J/2
    Output[2:(i-1)] = Delta_J
    Output[i]       = Delta_J/2 + (x-Knots[i])*(1+Knots[i]/Delta_J) - 0.5*(x^2-Knots[i]^2)/Delta_J
    Output[i+1]     = (x-Knots[i])*(1-Knots[i+1]/Delta_J) + 0.5*(x^2-Knots[i]^2)/Delta_J
  }
  return(Output)
}
