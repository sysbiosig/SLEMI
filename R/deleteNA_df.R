#' Removing NAs observations from data frame
#' 
#' @param data is a data.frame object
#' @return a data.frame object with the same structure as data and no observation with missing (NA) values
#' @export
#' @examples 
#' 
#' @keywords internal
deleteNA_df<-function(data){
  data[apply(data,1,function(x) !any(is.na(x)) ),]
  }