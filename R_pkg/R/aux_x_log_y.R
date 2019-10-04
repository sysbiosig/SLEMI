#' Calculation of expression \eqn{x\cdot \log y}
#' 
#' Internal, auxillary functions
#'
#' @param x is a numeric vector
#' @param y is a numeric vector (the same length of x)
#' @return Function calculates the value of expression \eqn{x\cdot \log y} element-wise in a numerically stable way. 
#' The result is a numeric vector of the same length as x. It is assumed that \eqn{0\cdot \log 0 = 0}.
#' @examples 
#' SLEMI:::aux_x_log_y(1,2)
#' SLEMI:::aux_x_log_y(0,0)
#' SLEMI:::aux_x_log_y(1000,100)
#' @keywords internal
#'
aux_x_log_y<-function(x,y){
  out=log(y^x)
  ids=is.infinite(out)
  out[ids]=x[ids]*log(y[ids])

  if (length(x)>100000){
  	temp=e1071::moment(x)
  	temp_str=Hmisc::capitalize("large")
  }

  out
}