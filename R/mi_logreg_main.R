#' Calculate mutual information
#' 
#' The main wrapping function for basic usage of SLEMI package for estimation of mutual information
#' @param dataRaw must be a data.frame object
#' @param pinput is a numeric vector with the probabilities of the input that will be used in estimation
#' @param signal is a character object that indicates columns of data to be treated as channel's input
#' @param response is a character vector that indicates columns of data to be treated as channel's output
#' @param side_variables is an optional character vector that indicates side variables' columns of data, if NULL no side variables are included
#' @param formula_string is a character object that includes a formula syntax to use in algorithm for capacity calculation
#' @param lr_maxit is the number of iteration of iterative algorithm of logistic regression
#' @param maxNWts is the maximum acceptable number of weights in logistic regression algorithm - the higher, the slower the computations
#' @param output_path is the directory in which output will be saved
#' @param plot_height, plot_width - the basic dimnesions of plots
#' @param model_out is the logical indicating if the calculated logisitc regression model should be included in output list
#' @param scale is the logical indicating if the data preprocessing should be carried out
#' @param graphs is the logical indicating if output graphs should be created
#' @param testing is the logical indicating if the testing procedures should be included
#' @param TestingSeed is the seed used in bootstrap procedures (only if testing=TRUE)
#' @param boot_num is the number of bootstrap tests to be performed (only if testing=TRUE)
#' @param boot_prob is the proportion of initial size of data to be used in bootstrap (only if testing=TRUE)
#' @param testing_cores - number of cores to use in parallel computing if implemented (only if testing=TRUE)
#' @param resamp_num is the number of resmapling tests to be performed (only if testing=TRUE)
#' @param traintest_num is the number of traintest tests to be performed (only if testing=TRUE)
#' @param partition_trainfrac is the fraction of data to be used as a training dataset (only if testing=TRUE)
#' @param glmnet_algorithm is the logical indicating if the glmnet package should be used
#' @param dataMatrix is a numeric matrix with columns treated as explanatory variables in 
#' logistic regression algorithm (only if glmnet_algorithm=TRUE)
#' @param glmnet_cores is the number of cores to use in parallel computing of glmnet package 
#' @param glmnet_lambdanum is the lambda parameter of glmnet package
#'
#' @return a list with several elements:
#' \itemize{
#' \item output$regression - confusion matrix of logistic regression predictions
#' \item output$cc         - channel capacity in bits
#' \item output$model      - nnet object describing logistic regression model (if model_out=TRUE)
#' \item output$params     - parameters used in algorithm
#' \item output$time       - computation time of calculations
#' \item output$testing    - a 4-element output list of testing procedures (bootstrap, resamplingMorph, traintest, bootResampMorph) (if testing=TRUE)
#' \item output$testing_pv - one-sided p-values of testing procedures (if testing=TRUE)
#' \item output$data       - dataRaw object
#' \item output$logGraphs  - ggplot and grid objects of output graphs (if graphs=TRUE)
#' }
#' @export
#' @examples 
#' tempdata=data_example1
#' dir.create("example1/",recursive=TRUE)
#' outputCLR1=mi_logreg_main(dataRaw=tempdata,
#' signal="signal", response="output",
#' lr_maxit=1500, output_path="example1/",plot_height=8,plot_width=12)
#' 
#' 
#' tempdata=data_example1
#' dir.create("example1_testing/",recursive=TRUE)
#' outputCLR1_testing=mi_logreg_main(dataRaw=tempdata,
#' signal="signal", response="output",
#' lr_maxit=1500, output_path="example1_testing/",plot_height=8,plot_width=12,
#' testing=TRUE,graphs=TRUE,TestingSeed=11111, boot_num=50,boot_prob=0.8,testing_cores=2,
#' sidevar_num=2,traintest_num=50,partition_trainfrac=0.6)
#' 
#' 
#' tempdata=data_example2
#' dir.create("example2/",recursive=TRUE)
#' outputCLR2=mi_logreg_main(dataRaw=tempdata,
#' signal="signal", response=c("X1","X2","X3"),
#' lr_maxit=1500, output_path="example2/",plot_height=8,plot_width=12) 
#' 
#' For further details see vignette
mi_logreg_main<-function(dataRaw, signal="signal", response="response",side_variables=NULL,
                         pinput=NULL, 
                                          formula_string=NULL,
                                          glmnet_algorithm=FALSE,dataMatrix=NULL, 
                                          glmnet_cores=1,glmnet_lambdanum=10,
                                          lr_maxit=1000, MaxNWts = 5000, 
                                          output_path=NULL,
                                          testing=FALSE, model_out=TRUE,scale=TRUE,graphs=TRUE,
                                          TestingSeed=1234,testing_cores=1,
                                          boot_num=10,boot_prob=0.8,
                                          sidevar_num=10,
                                          traintest_num=10,partition_trainfrac=0.6,
                                          plot_width=6,plot_height=4,
                                          dataout=TRUE){
  
  #Debugging:
  print("Procedure starting")
  
  time_start=proc.time()
  
  # checking assumptions
  if (is.null(output_path)) { 
    warning('path is not defined. Graphs and RDS file will not be saved.')
    }
  if (!is.data.frame(dataRaw)) {
    stop('data is not in data.frame format')
  }
  if ( length(colnames(dataRaw)==signal)==0 ) {
    stop('There is no column described as signal in data')
  }
  if ( length(colnames(dataRaw) %in% response)==length(response) ) {
    stop('There is no column described as response in data')
  }
  if (!is.null(side_variables)){
    if ( length(colnames(dataRaw) %in% side_variables)==length(side_variables) ) {
    stop('There is no column described as side_variables in data')
    }
  }

  
   data0=dataRaw[,c(signal,response,side_variables)]
   if ( any(apply(data0,1,function(x) any(is.na(x)) )) ) {
     print("There are NA in observations - removing")
     data0=data0[!apply(data0,1,function(x) any(is.na(x)) ),]
     print("Numer of observations after cleaning:")
     print(table(data0[[signal]]))
   }
  
   data0=aux_signal_transform(data0,signal)
   tempcolnames=colnames(data0)
   tempsignal=data.frame(data0[,(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )])
   colnames(tempsignal)<-tempcolnames[(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )]
   data0=data.frame(data0[,!(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )])
   colnames(data0)<-tempcolnames[!(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )]
   
  print("Preprocessing started")
  
  #PreProcessing
  temp_idnumeric=sapply(data0,is.numeric)
  if (scale&sum(temp_idnumeric)==1) {
    data0[,temp_idnumeric]<-(data0[,temp_idnumeric]-mean(data0[,temp_idnumeric]))/sd(data0[,temp_idnumeric])
    data <- cbind(data0,tempsignal)
  } else if (scale) {
    preProcValues <- caret::preProcess(data0, method = c("center", "scale"))
    data <- cbind(predict(preProcValues, data0),tempsignal)
  } else {
    data <- cbind(data0,tempsignal)
  }
  rm(temp_idnumeric)

  #Debugging:
  print("Preprocessing completed. Algorithm initialization")
  

  
  output<-mi_logreg_algorithm(data=data,signal=signal,response=response,side_variables=side_variables,
                              pinput=pinput,
                                              formula_string=formula_string, model_out = model_out,
                                              lr_maxit=lr_maxit,MaxNWts =MaxNWts) 

  
  #Debugging:
  print(paste("Main algorihtm sent data:",length(output)))
  
  if (testing){
      output$testing<-mi_logreg_testing(data,signal=signal,response=response,side_variables=side_variables,
                                                      lr_maxit=lr_maxit,MaxNWts =MaxNWts,
                                                      formula_string=formula_string, model_out = FALSE,
                                                      TestingSeed=TestingSeed,testing_cores=testing_cores,
                                                      boot_num=boot_num,boot_prob=boot_prob,
                                                      sidevar_num=sidevar_num,
                                                      traintest_num=traintest_num,partition_trainfrac=partition_trainfrac,
                                              glmnet_algorithm= glmnet_algorithm,dataMatrix=dataMatrix, 
                                              glmnet_lambdanum=glmnet_lambdanum
                                              )
      
      output$testing_pv<-lapply(output$testing,function(x){
        tmp_boot_cc=sapply(x,function(xx) xx$cc)
        c(mean(tmp_boot_cc<output$cc),mean(tmp_boot_cc>output$cc))
      })
  }
  
  
  #Debugging:
  print(paste("Testing algorithm sent data:", length(output$testing)))
  

  output$time   <- proc.time() - time_start
  output$params <- c(lr_maxit=lr_maxit,MaxNWts =MaxNWts)
  
  if (dataout){
  output$data   <- dataRaw
  }
  
  if(!is.null(output_path)){dir.create(output_path,recursive=TRUE)}
  
  
  #Debugging:
  print("RDS saved")
  
  if (graphs){
    output$logGraphs=try(capacity_output_graphs(data=dataRaw,signal=signal,response=response,side_variables=side_variables,cc_output=output,
                                output_path=output_path,height=plot_height,width=plot_width),
        silent=FALSE)
    
    #Debugging:
    print("Graphs finished")
  }
  
  output$mi=output$cc
  output$cc=NULL
  
  if(!is.null(output_path)){saveRDS(output,file=paste(output_path,'output.rds',sep=""))}
  
  output
}