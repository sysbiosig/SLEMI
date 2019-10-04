#' Estimate channel capacity between discrete input and continuouse output
#'
#' The main wrapping function for basic usage of SLEMI package for estimation of channel capacity. Firstly, data is pre-processed
#' (all arguments are checked, observation with NAs are removed, variables are scaled and centered (if scale=TRUE)). Then basic estimation is carried out
#' and (if testing=TRUE) diagnostic tests are computed. If output directory path is given (output_path is not NULL), graphs visualising the data and the analysis
#' are saved there, together with a compressed output object (as .rds file) with full estimation results.
#'
#' In a typical experiment aimed to quantify information flow a given signaling system, input values \eqn{x_1\leq x_2 \ldots... \leq x_m}, ranging from 0 to saturation are considered.
#' Then, for each input level, \eqn{x_i}, \eqn{n_i} observations are collected, which are represetned as vectors 
#' \deqn{y^i_j \sim P(Y|X = x_i)}
#' Within information theory the degree of information transmission is measured as the mutual information
#' \deqn{MI(X,Y) = \sum_{i=1}^{m} P(x_i)\int_{R^k} P(y|X = x_i)log_2\frac{P(y|X = x_i)}{P(y)}dy,}
#' where \eqn{P(y)} is the marginal distribution of the output. MI is expressed in bits and \eqn{2^{MI}} can be interpreted as the number of 
#' inputs that the system can resolve on average.
#'
#' The maximization of mutual information with respect to the input distribution, \eqn{P(X)},
#' defines the information capacity, C. Formally,
#' \deqn{C^* = max_{P(X)} MI(X,Y)}
#' Information capacity is expressed in bits and \eqn{2^{C^*}} can be interpreted as the maximal number of inputs that the system can
#' effectively resolve.
#'
#' In contrast to existing approaches, instead of estimating, possibly highly dimensional, conditional output distributions P(Y|X =x_i), we propose to estimate the discrete, conditional input distribution, 
#' \eqn{P(x_i |Y = y)}, which is known to be a simpler problem. Estimation of the MI using estimates of \eqn{P(x_i |Y = y)}, denoted here as \eqn{\hat{P}(x_i|Y = y)}, is possible as the MI, can be
#' alternatively written as
#' \deqn{MI(X,Y) = \sum_{i=1}^{m} P(x_i)\int_{R^k} P(y|X = x_i)log_2\frac{P(x_i|Y = y)}{P(x_i)}dy}
#' The expected value (as in above expression) with respect to distribution \eqn{P(Y|X = x_i)} can be approximated by the average with respect to data
#' \deqn{MI(X,Y) \approx \sum_{i=1}^{m} P(x_i)\frac{1}{n_i} \sum_{j=1}^{n_i} P(y|X = x_i)log_2\frac{\hat{P}(x_i|Y = y^i_j)}{P(x_i)}dy}
#' Here, we propose to use logistic regression as \eqn{\hat{P}(x_i|Y = y^i_j)}. Specifically,
#' \deqn{log\frac{P(x_i |Y = y)}{P(x_m|Y = y)} \approx \alpha_i +\beta_iy}
#'
#' Following this approach, channel capacity can be calculated by optimising MI with respect to the input distribution, \eqn{P(X)}.
#' However, this, potentially difficult problem, can be divided into two simpler maximization problems, for which explicit solutions exist. 
#' Therefore, channel capacity can be obtained from the two explicit solutions in an iterative procedure known as alternate maximization (similarly as in Blahut-Arimoto algorithm) [1].
#'
#' Additional parameters: lr_maxit and maxNWts are the same as in definition of multinom function from nnet package. An alternative
#' model formula (using formula_string arguments) should be provided if  data are not suitable for description by logistic regression
#' (recommended only for advanced users). Preliminary scaling of  data (argument scale) should be used similarly as in other 
#' data-driven approaches, e.g. if response variables are comparable, scaling (scale=FALSE) can be omitted, while if they represent 
#' different phenomenon (varying by units and/or magnitude) scaling is recommended.
#'
#' @section References:
#' [1] Csiszar I, Tusnady G, Information geometry and alternating minimization procedures, Statistics & Decisions 1 Supplement 1 (1984), 205–237.
#'
#' [2] Jetka T, Nienaltowski K, Winarski T, Blonski S, Komorowski M,  
#' Information-theoretic analysis of multivariate single-cell signaling responses using SLEMI,
#' \emph{PLoS Comput Biol}, 15(7): e1007132, 2019, https://doi.org/10.1371/journal.pcbi.1007132.
#'
#'
#'
#' @param dataRaw must be a data.frame object
#' @param signal is a character object with names of columns of dataRaw to be treated as channel's input.
#' @param response is a character vector with names of columns of dataRaw  to be treated as channel's output
#' @param output_path is the directory in which output will be saved
#' @param side_variables (optional) is a character vector that indicates side variables' columns of data, if NULL no side variables are included
#' @param formula_string (optional) is a character object that includes a formula syntax to use in logistic regression model. 
#' If NULL, a standard additive model of response variables is assumed. Only for advanced users.
#' @param cc_maxit is the number of iteration of iterative optimisation of the algorithm to esimate channel capacity. Default is 100.
#' @param lr_maxit is a maximum number of iteration of fitting algorithm of logistic regression. Default is 1000.
#' @param MaxNWts is a maximum acceptable number of weights in logistic regression algorithm. Default is 5000.
#' @param plot_height - the basic dimnesions (height) of plots, in inches
#' @param plot_width - the basic dimnesions (width) of plots, in inches
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
#' outputCLR1=capacity_logreg_main(dataRaw=tempdata,
#' signal="signal", response="response",
#' formula_string = "signal~response")
#' 
#' tempdata=data_example2
#' outputCLR2=capacity_logreg_main(dataRaw=tempdata,
#' signal="signal", response=c("X1","X2","X3"),
#' formula_string = "signal~X1+X2+X3") 
#' 
#' #For further details see vignette
capacity_logreg_main<-function(dataRaw, signal="input", response=NULL,output_path=NULL,
  side_variables=NULL,formula_string=NULL,cc_maxit=100,lr_maxit=1000, MaxNWts = 5000, 
  testing=FALSE, model_out=TRUE,scale=TRUE,
                                          TestingSeed=1234,testing_cores=1,
                                          boot_num=10,boot_prob=0.8,
                                          sidevar_num=10,
                                          traintest_num=10,partition_trainfrac=0.6,
                                          plot_width=6,plot_height=4,
                                          data_out=FALSE){
  
 if(!is.null(output_path)){
    dir.create(output_path,recursive=TRUE)
 }
  
  #Debugging:
  message(" Estimating channel capacity ...")
  
   if (is.null(response)){
    response=paste0("output_",1:(ncol(dataRaw)-1) )  
  }
  
  time_start=proc.time()
  dataRaw=as.data.frame(dataRaw)

  # checking assumptions
  if (is.null(output_path)) { 
    message('Output path is not defined. Graphs and RDS file will not be saved.')
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
     message("There are NA in observations - removing...")
     data0=data0[!apply(data0,1,function(x) any(is.na(x)) ),]
   }
  
   data0=func_signal_transform(data0,signal)
   tempcolnames=colnames(data0)
   tempsignal=data.frame(data0[,(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )])
   colnames(tempsignal)<-tempcolnames[(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )]
   data0=data.frame(data0[,!(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )])
   colnames(data0)<-tempcolnames[!(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )]
   
 # message(" Preprocessing started")
  
  #PreProcessing
  temp_idnumeric=sapply(data0,is.numeric)
  if (scale&sum(temp_idnumeric)==1) {
    data0[,temp_idnumeric]<-(data0[,temp_idnumeric]-mean(data0[,temp_idnumeric]))/stats::sd(data0[,temp_idnumeric])
    data <- cbind(data0,tempsignal)
  } else if (scale) {
    preProcValues <- caret::preProcess(data0, method = c("center", "scale"))
    data <- cbind(stats::predict(preProcValues, data0),tempsignal)
  } else {
    data <- cbind(data0,tempsignal)
  }
  rm(temp_idnumeric)

  #Debugging:
  #message("Running main algorithm...")
  
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
    dir.create(output_path,recursive=TRUE)

    message(" Drawing graphs and saving objects ...")
    temp_logGraphs=try(output_graphs_main(data=dataRaw,signal=signal,response=response,side_variables=side_variables,cc_output=output,
                                output_path=output_path,height=plot_height,width=plot_width),
        silent=FALSE)
    rm(temp_logGraphs)
    #Debugging:

    saveRDS(output,file=paste(output_path,'output.rds',sep=""))
    message(paste0(" Estimation finished. Results saved in ",output_path,""))  
  } else {
    message(paste0(" Estimation finished."))
  }
  
  output
}
