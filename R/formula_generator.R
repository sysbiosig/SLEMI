#' Formula for logistic regression algorithm
#' 
#' @param signal is a character object that indicates columns of data to be treated as channel's input
#' @param response is a character vector that indicates columns of data to be treated as channel's output
#' @param side_variables is a character vector that indicates side variables' columns of data
#' @return A character object that includes a standard formula syntax to use in algorithm for capacity calculation
#' @export
#' @examples 
#' formula_generator(signal="signal",response="response", side_variables=NULL)
#' formula_generator(signal="inputX",response="responseY", side_variables="SV1")
#' formula_generator(signal="signalX",response=c("r_1","r_2","r_5"), side_variables="SV")
formula_generator<-function(signal="signal",response="response", side_variables=NULL){
  if (!is.null(side_variables)){
    formula_string<- (paste(signal,"~",
                            paste(
                              paste(c(response),collapse="+"),
                              paste(apply(expand.grid(side_variables,response),1,function(x) paste(x,collapse=":")),collapse="+"),
                              sep="+"),
                            collapse=""))
  } else {
    formula_string<- (paste(signal,"~",
                            paste(
                              paste(c(response),collapse="+"),
                              sep="+"),
                            collapse=""))
  }
  formula_string
}
