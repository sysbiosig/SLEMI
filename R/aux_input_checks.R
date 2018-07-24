#' Auxiliary function. Initial verification of input
#' 
#' @param data is an input object that should be a data.frame
#' @param signal is a character object that indicates input columns of data
#' @param response is a character vector that indicates output's columns of data
#' @param side_variables is a character vector that indicates side variables' columns of data
#' @return If all initial data is valid, string "ok" is retured. Otherwise, error is given.
#' @keywords internal
#' @export
#' @examples 
#' data=data_example1
#' aux_input_checks(data=data,signal="signal",response="response",side_variables="sideVar")
#' 
#' data=as.matrix(data_example1)
#' aux_input_checks(data=data,signal="signal",response="response",side_variables="sideVar")
#' 
#' data=data_example1
#' aux_input_checks(data=data,signal="input",response="response",side_variables="sideVar")
aux_input_checks<-function(data,signal,response,side_variables){
  if (!is.data.frame(data)) {
    stop('data is not in data.frame format')
  }
  if ( length(colnames(data)==signal)==0 ) {
    stop('There is no column described as signal in data')
  }
  if ( length(colnames(data)%in%response)==length(response) ) {
    stop('There is no column described as response in data')
  }
  if ( length(colnames(data)%in%side_variables)==length(side_variables) ) {
    stop('There is no column described as side_variables in data')
  }
  if ( any(apply(data,1,function(x) any(is.na(x)) )) ) {
    stop("There are NA in observations")
  }
  out="ok"
}