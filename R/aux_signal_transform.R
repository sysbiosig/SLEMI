#' Auxiliary function. Initial verification of signal provided and its transformation into a factor type
#' 
#' @param data is a data.frame
#' @param signal is a character that indicates columns of data that include the labels of input
#' @return A data.frame that is a copy of data provided with signal column transformed to factor class. 
#' If signal has been numeric initially, additional column is created "signal_RAW" that is an exact
#' copy of original column
#' @export
#' @examples 
#' data=data_example1
#' data1=aux_signal_transform(data,"signal")
#' data$signal=as.character(data$signal)
#' data2=aux_signal_transform(data,"signal")
#' data$signal=as.numeric(data$signal)
#' data3=aux_signal_transform(data,"signal")
#' @keywords internal
aux_signal_transform<-function(data,signal){
  signal_class=class(data[[signal]])
  
    if (signal_class=="numeric") {
      data[[paste(signal,"_RAW",sep="")]]=data[[signal]]
      data[[signal]]=factor(data[[signal]],levels=sort(unique(data[[signal]])))
    } else if (signal_class=="character") {
      tmp_signal_values=unique(data[[signal]])
      tmp_signal_values_num=as.numeric(tmp_signal_values)
    if (any(is.na(tmp_signal_values_num))){
      tmp_signal_values_num=stringr::str_match(tmp_signal_values,"[0-9.]+")
      if (any(is.na(tmp_signal_values_num))) {
        data[[signal]]=factor(data[[signal]])
      } else {
        data[[signal]]=as.numeric(stringr::str_match(data[[signal]],"[0-9.]+"))
        data[[paste(signal,"_RAW",sep="")]]=data[[signal]]
        data[[signal]]=factor(data[[signal]],levels=sort(unique(data[[signal]])))
      }
    } else {
      data[[signal]]=factor(as.numeric(data[[signal]]),levels=sort(unique(as.numeric(data[[signal]]) ))) 
    } 
  } else if (signal_class=="factor"){
    print("")
  } else {
    stop('signal is of unknown type')
  }
    
  data
}