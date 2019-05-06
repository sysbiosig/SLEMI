#' Estimate channel capacity between discrete input and continuouse output
#'
#' The main wrapping function for basic usage of SLEMI package for estimation of channel capacity. Firstly, data is pre-processed
#' (all arguments are checked, observation with NAs are removed, variables are scaled and centered (if scale=TRUE)). Then basic estimation is carried out
#' and (if testing=TRUE) diagnostic tests are computed. If output directory path is given (output_path is not NULL), graphs visualising the data and the analysis
#' are saved there, together with a compressed output object (as .rds file) with full estimation results.
#'
#' Additional parameters: lr_maxit and maxNWts are the same as in definition of multinom function from nnet package. An alternative
#' model formula (using formula_string arguments) should be provided if  data are not suitable for description by logistic regression
#' (recommended only for advanced users). Preliminary scaling of  data (argument scale) should be used similarly as in other 
#' data-driven approaches, e.g. if response variables are comparable, scaling (scale=FALSE) can be omitted, while if they represent 
#' different phenomenon (varying by units and/or magnitude) scaling is recommended.
#'
#' @section References:
#' Jetka T, Nienaltowski K, Winarski T, Blonski S, Komorowski M,  
#' Information-theoretic analysis of multivariate single-cell signaling responses using SLEMI,
#' \emph{PLOS Comp Bio}, 2019.
#'
#' @param dataRaw must be a data.frame object
#' @param signal is a character object with names of columns of dataRaw to be treated as channel's input.
#' @param response is a character vector with names of columns of dataRaw  to be treated as channel's output
#' @param side_variables (optional) is a character vector that indicates side variables' columns of data, if NULL no side variables are included
#' @param formula_string (optional) is a character object that includes a formula syntax to use in logistic regression model. 
#' If NULL, a standard additive model of response variables is assumed. Only for advanced users.
#' @param cc_maxit is the number of iteration of iterative optimisation of the algorithm to esimate channel capacity. Default is 100.
#' @param lr_maxit is a maximum number of iteration of fitting algorithm of logistic regression. Default is 1000.
#' @param maxNWts is a maximum acceptable number of weights in logistic regression algorithm. Default is 5000.
#' @param output_path is the directory in which output will be saved
#' @param plot_height, plot_width - the basic dimnesions of plots
#' @param model_out is the logical indicating if the calculated logisitc regression model should be included in output list
#' @param data_out  is the logical indicating if the data should be included in output list
#' @param scale is a logical indicating if the response variables should be scaled and centered before fitting logistic regression
#' @param testing is the logical indicating if the testing procedures should be executed
#' @param TestingSeed is the seed for random number generator used in testing procedures
#' @param boot_num is the number of bootstrap tests to be performed. Default is 10, but it is recommended to use at least 50 for reliable estimates.
#' @param boot_prob is the proportion of initial size of data to be used in bootstrap
#' @param testing_cores - number of cores to be used in parallel computing (via doParallel package)
#' @param sidevar_num is the number of re-shuffling tests of side variables to be performed. Default is 10, but it is recommended to use at least 50 for reliable estimates.
#' @param traintest_num is the number of overfitting tests to be performed. Default is 10, but it is recommended to use at least 50 for reliable estimates.
#' @param partition_trainfrac is the fraction of data to be used as a training dataset
#'
#' @return a list with several elements:
#' \itemize{
#' \item output$regression - confusion matrix of logistic regression predictions
#' \item output$cc         - channel capacity in bits
#' \item output$p_opt      - optimal probability distribution
#' \item output$model      - nnet object describing logistic regression model (if model_out=TRUE)
#' \item output$params     - parameters used in algorithm
#' \item output$time       - computation time of calculations
#' \item output$testing    - a 2- or 4-element output list of testing procedures (if testing=TRUE)
#' \item output$testing_pv - one-sided p-values of testing procedures (if testing=TRUE)
#' \item output$data       - raw data used in analysis
#' }
#' @export
#' @examples 
#' tempdata=data_example1
#' dir.create("example1/",recursive=TRUE)
#' outputCLR1=capacity_logreg_main(dataRaw=tempdata,
#' signal="signal", response="response",
#' formula_string = "signal~output",
#' cc_maxit=75,lr_maxit=1500, output_path="example1/",plot_height=8,plot_width=12)
#' 
#' 
#' tempdata=data_example1
#' dir.create("example1_testing/",recursive=TRUE)
#' outputCLR1_testing=capacity_logreg_main(dataRaw=tempdata,
#' signal="signal", response="response",
#' cc_maxit=75,lr_maxit=1500, output_path="example1_testing/",plot_height=8,plot_width=12,
#' testing=TRUE,graphs=TRUE,TestingSeed=11111, boot_num=50,boot_prob=0.8,testing_cores=2,
#' sidevar_num=2,traintest_num=50,partition_trainfrac=0.6)
#' 
#' 
#' tempdata=data_example2
#' dir.create("example2/",recursive=TRUE)
#' outputCLR2=capacity_logreg_main(dataRaw=tempdata,
#' signal="signal", response=c("X1","X2","X3"),
#' formula_string = "signal~X1+X2+X3",
#' cc_maxit=75,lr_maxit=1500, output_path="example2/",plot_height=8,plot_width=12) 
#' 
#' For further details see vignette
capacity_logreg_main<-function(dataRaw, signal="input", response=NULL,side_variables=NULL,
                                          formula_string=NULL,
                                          cc_maxit=100,lr_maxit=1000, MaxNWts = 5000, 
                                          output_path=NULL,
                                          testing=FALSE, model_out=TRUE,scale=TRUE,
                                          TestingSeed=1234,testing_cores=1,
                                          boot_num=10,boot_prob=0.8,
                                          sidevar_num=10,
                                          traintest_num=10,partition_trainfrac=0.6,
                                          plot_width=6,plot_height=4,
                                          data_out=FALSE){
  
  #Debugging:
  cat("\n Estimating channel capacity ...")
  
   if (is.null(response)){
    response=paste0("output_",1:(ncol(dataRaw)-1) )  
  }
  
  time_start=proc.time()
  dataRaw=as.data.frame(dataRaw)

  # checking assumptions
  if (is.null(output_path)) { 
    warning('path is not defined. Graphs and RDS file will not be saved.')
    }
  if (!is.data.frame(dataRaw)) {
    stop('data is not in data.frame format')
  }
  if ( sum(colnames(dataRaw)==signal)==0 ) {
    stop('There is no column described as signal in data')
  }
  if (!sum(colnames(dataRaw) %in% response)==length(response) ) {
    stop('There is no column described as response in data')
  }
  if (!is.null(side_variables)){
    if (!sum(colnames(dataRaw) %in% side_variables)==length(side_variables) ) {
    stop('There is no column described as side_variables in data')
    }
  }

  
   data0=dataRaw[,c(signal,response,side_variables)]
   if ( any(apply(data0,1,function(x) any(is.na(x)) )) ) {
     cat("\nThere are NA in observations - removing...")
     data0=data0[!apply(data0,1,function(x) any(is.na(x)) ),]
     cat("Numer of observations after cleaning:")
     cat(table(data0[[signal]]))
   }
  
   data0=func_signal_transform(data0,signal)
   tempcolnames=colnames(data0)
   tempsignal=data.frame(data0[,(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )])
   colnames(tempsignal)<-tempcolnames[(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )]
   data0=data.frame(data0[,!(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )])
   colnames(data0)<-tempcolnames[!(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )]
   
  cat("\n Preprocessing started")
  
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
  cat("... Preprocessing completed.")
  
  output<-capacity_logreg_algorithm(data=data,signal=signal,response=response,side_variables=side_variables,
                                              formula_string=formula_string, model_out = model_out,
                                              cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts) 
  
  #Debugging:
  #print(paste("Main algorihtm sent data:",length(output)))
  
  if (testing){
      output$testing<-capacity_logreg_testing(data,signal=signal,response=response,side_variables=side_variables,
                                                      cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts,
                                                      formula_string=formula_string, 
                                                      TestingSeed=TestingSeed,testing_cores=testing_cores,
                                                      boot_num=boot_num,boot_prob=boot_prob,
                                                      sidevar_num=sidevar_num,
                                                      traintest_num=traintest_num,partition_trainfrac=partition_trainfrac)

      output$testing_pv<-lapply(output$testing,function(x){
        tmp_boot_cc=sapply(x,function(xx) xx$cc)
        c(mean(tmp_boot_cc<output$cc),mean(tmp_boot_cc>output$cc))
      })
  }
  
  
  #Debugging:
  #print(paste("Testing algorithm sent data:", length(output$testing)))
  
  output$time   <- proc.time() - time_start
  output$params <- c(cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts)
  
  if (data_out){
    output$data   <- dataRaw
  }
  
  if(!is.null(output_path)){
    options(warn=-1)
    dir.create(output_path,recursive=TRUE)
    options(warn=0)

    cat("\n Creating graphs ...")
    temp_logGraphs=try(output_graphs_main(data=dataRaw,signal=signal,response=response,side_variables=side_variables,cc_output=output,
                                output_path=output_path,height=plot_height,width=plot_width),
        silent=FALSE)
    rm(temp_logGraphs)
    #Debugging:
    cat("finished")

    saveRDS(output,file=paste(output_path,'output.rds',sep=""))
    cat(paste0("\n Full Procedure finished. Results are in ",output_path))  
  } else {
    cat(paste0("\n Full Procedure finished."))
  }
  
  output
}
