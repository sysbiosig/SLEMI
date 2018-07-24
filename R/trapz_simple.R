#' Numerical intergration by trapezoidal rule
#' 
#' This function calculate the integral of function described by set of points (x,y) using trapezoidal rule.
#' @param x is a numeric vector
#' @param y is a numeric vector of the same length as x
#' @return A numeric value of integral approximated by trapezoids
#' @export
#' @examples 
#' trapz_simple(x=1:10,y=exp(1:10))
#' @keywords internal
trapz_simple<-function(x,y){
  nx=length(x)
  
  xmesh=x[(2:nx)]-x[(1:(nx-1))]
  ymesh=0.5*(y[(2:nx)]+y[(1:(nx-1))])
  
  sum(xmesh*ymesh)
}