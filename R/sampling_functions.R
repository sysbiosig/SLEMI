#' Sampling procedures used for testing capacity algorithm
#' 
#' These function allow to re-sample, bootstrap and divide initial dataset
#' @param data is a data.frame to be resampled
#' @param dataDiv a character indicating column of data, with respect to which, data should be splitte before bootstrap
#' @param prob is numeric for the portion of data that should be sampled from the whole dataset (only in sampling_bootstrap)
#' @param side_variables is a vector of characters indicating columns of data the will be reshuffled (only in sampling_shuffle)
#' @param partition_trainfrac is a numeric for the portion of data that will be used as a training and testing datasets
#' 
#' @return Function sampling_bootstrap returns a data.frame with the same structure as initial data object, but with prob proportion
#' of observations for each dataDiv level. Function sampling_shuffle returns a data.frame with the same structure as initial data object with 
#' shuffled values of columns given in side_variables argument. Function sampling_partition returns a list of two data.frame objects - 
#' train and test that has the same structure as initial data argument with partition_trainfrac and 1-partition_trainfrac observations, resepctively.

#' @export
#' @examples 
#' data=data_example1
#' dataBootstrap = sampling_bootstrap(data=data,prob=0.8,"signal")
#' dataShuffle = sampling_shuffle(data=data,"sideVar")
#' dataTrainTest = sampling_partition(data=data,dataDiv="signal",partition_trainfrac=0.6)

sampling_bootstrap<-function(data,prob=1,dataDiv){
  # Bootstrapping samples within a given variable level
  # The output data has ~ prob*ObsNum observations
  #   signal="signal"
  #   colnames(data)[colnames(data)==var]="signal"
  
  dataBoot=do.call(rbind,by(data,list(dataDiv),function(x){ 
    sample_num=nrow(x)
    sample_out=sample(1:sample_num,floor(prob*sample_num),replace=FALSE)
    x[sample_out,]
  }))
  
  #   levels_data=levels(data[,signal])
  #   datanew=data
  #   for (k in 1:length(levels_data)){
  #     sample_num=sum(datanew[,var]==levels_data[k])
  #     sample_out=sample(1:sample_num,floor(prob*sample_num),replace=TRUE)
  #     IDout<-(1:nrow(datanew))[datanew[,var]==levels_data[k]][sample_out]
  #     datanew=datanew[-IDout,]
  #   }
  #   datanew
  colnames(dataBoot)<-colnames(data)
  if (class(data)=="matrix") {dataBoot=as.matrix(dataBoot)}
  
  dataBoot
}


#' @rdname sampling_bootstrap
sampling_shuffle<-function(data,side_variables){
  
  if (!is.null(side_variables)){
    
    if (is.matrix(data)) {
      dataBasic=data[,side_variables]
    } else if (!is.data.frame(data)){
      stop('wrong data format')
  #  } else if (is.data.table(data)) {
  #   dataBasic=data[,!(colnames(data)%in%side_variables),with=FALSE]
    } else {dataBasic=data[,!(colnames(data)%in%side_variables)]}
    
    obsv_num=nrow(data)
    
    idsamples=replicate(length(side_variables),sample(1:obsv_num,obsv_num))
    dataNewAdd=do.call(cbind,lapply(1:length(side_variables),function(k){
      if (is.matrix(data))  { 
        datax=data[idsamples[,k],side_variables[k]]
      } else if (!is.data.frame(data)){
        stop('wrong data format')
   #   } else if (is.data.table(data)) {
   #     datax=data[idsamples[,k],side_variables[k],with=FALSE]
      } else {datax=data[idsamples[,k],side_variables[k]]}
    }))
    
    colnames(dataNewAdd)=side_variables
    dataNew=cbind(dataBasic,dataNewAdd)
    if (class(data)=="matrix") {dataNew=as.matrix(dataNew)}
    
  } else {
    warning('there is no side variables defined. No resampling made.')
    dataNew=data
  }
  
  dataNew
}


#' @rdname sampling_bootstrap
sampling_partition<-function(data,dataDiv,partition_trainfrac){
  dataNew=list()
  
  #   colnames(data)[colnames(data)==signal]="signal"
  #   signal="signal"
  
  group_partition=by(data,list(dataDiv),function(x){ 
    outx=list()
    sample_num=nrow(x)
    idsample<-sample(1:sample_num,floor(partition_trainfrac*sample_num))
    outx$datatrain=x[idsample,]
    outx$datatest=x[-idsample,]
    outx
  })
  
  dataNew$test=do.call(rbind,lapply(group_partition,function(x) x$datatest ))
  dataNew$train=do.call(rbind,lapply(group_partition,function(x) x$datatrain ))
  
  if (class(data)=="matrix") {
    dataNew$test=as.matrix(dataNew$test)
    dataNew$train=as.matrix(dataNew$train)}
  
  dataNew
}