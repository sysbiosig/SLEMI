#' Removing NAs observations from data frame
#' 
#' @param data is a data.frame object
#' @return a data.frame object with the same structure as data and no observation with missing (NA) values
#' @examples 
#' df=data.frame(x=c(rnorm(10),NA,NA),y=c(NA,NA,rnorm(10)))
#' aux_deleteNA_df(df)
#' @keywords internal
aux_deleteNA_df<-function(data){
  data[apply(data,1,function(x) !any(is.na(x)) ),]
}