#' Testing procedures for calculating mutual information
#' 
#' INPUT:
#' @param data must be a data.frame object
#' @param signal is a character object that indicates columns of data to be treated as channel's input
#' @param pinput is a numeric vector with the probabilities of the input that will be used in estimation
#' @param response is a character vector that indicates columns of data to be treated as channel's output
#' @param side_variables is an optional character vector that indicates side variables' columns of data, if NULL no side variables are included
#' @param formula_string is a character object that includes a formula syntax to use in algorithm for capacity calculation
#' @param lr_maxit is the number of iteration of iterative algorithm of logistic regression
#' @param maxNWts is the maximum acceptable number of weights in logistic regression algorithm - the higher, the slower the computations
#' @param TestingSeed is the seed used in bootstrap procedures (only if testing=TRUE)
#' @param boot_num is the number of bootstrap tests to be performed (only if testing=TRUE)
#' @param boot_prob is the proportion of initial size of data to be used in bootstrap (only if testing=TRUE)
#' @param testing_cores - number of cores to use in parallel computing if implemented (only if testing=TRUE)
#' @param resamp_num is the number of resmapling tests to be performed (only if testing=TRUE)
#' @param traintest_num is the number of traintest tests to be performed (only if testing=TRUE)
#' @param partition_trainfrac is the fraction of data to be used as a training dataset (only if testing=TRUE)
#' @keywords internal
#' @return a list with four elements:
#' \itemize{
#' \item output$bootstrap - confusion matrix of logistic regression predictions
#' \item output$resamplingMorph          - channel capacity in bits
#' \item output$traintest      - optimal probability distribution
#' \item output$bootResampMorph      - nnet object describing logistic regression model (if model_out=TRUE)
#' }
#' Each of above is a list, where an element is an output of a single repetition of the channel capacity algorithm
#' @export
#' @examples 
#' 
mi_logreg_testing<-function(data,signal="signal",response="response",side_variables=NULL,
                                  lr_maxit=1000,MaxNWts = 5000,
                                  formula_string=NULL,
                                  model_out=FALSE,
                                  TestingSeed=1234,testing_cores=4,
                                  boot_num=10,boot_prob=0.8,
                                  sidevar_num=10,
                                  traintest_num=10,partition_trainfrac=0.6,
                                  pinput=NULL){
  
  output=list()
  
  data_signal=data[[signal]]
  
  set.seed(TestingSeed)
  cat("Testing started..")
  
  `%dopar%`<-foreach::`%dopar%`
  
  #Bootstrap
  cl=parallel::makeCluster(testing_cores)
  doParallel::registerDoParallel(cl)
  output_test1<-foreach::foreach(j=1:boot_num,
                                 .export=c("sampling_bootstrap","mi_logreg_algorithm","aux_x_log_y"),
                                 .packages=c("nnet","caret")) %dopar% {
                                     data_bt_samp   <- sampling_bootstrap(data,boot_prob,data_signal)
                                     bt_samp_output <- mi_logreg_algorithm(data=data_bt_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                 formula_string=formula_string,model_out=model_out,lr_maxit=lr_maxit,MaxNWts =MaxNWts,
                                                                           pinput=pinput)
                                   bt_samp_output
                                 }
  parallel::stopCluster(cl)
  cat("Bootstrap completed..")
  
  if (!is.null(side_variables)){
    # Resampling
    cl=parallel::makeCluster(testing_cores)
    doParallel::registerDoParallel(cl)
    output_test2=foreach::foreach(j=1:sidevar_num,
                                  .export=c("sampling_shuffle","mi_logreg_algorithm","aux_x_log_y"),
                                  .packages=c("nnet","caret")) %dopar% {
                                      data_resamp_samp   <- sampling_shuffle(data,side_variables)
                                      bt_samp_output <- mi_logreg_algorithm(data=data_resamp_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                  formula_string=formula_string,model_out=model_out,lr_maxit=lr_maxit,MaxNWts =MaxNWts,
                                                                            pinput=pinput )
                                    bt_samp_output
                                  }
    parallel::stopCluster(cl)
    cat("Resampling completed..")
    
    
    # Bootstrap&Resampling
    cl=parallel::makeCluster(testing_cores)
    doParallel::registerDoParallel(cl)
    output_test4=foreach::foreach(j=1:sidevar_num,
                                  .export=c("sampling_bootstrap","sampling_shuffle","mi_logreg_algorithm","aux_x_log_y"),
                                  .packages=c("nnet","caret")) %dopar% {
                                      data_resamp_samp   <- sampling_shuffle(data,side_variables)
                                      data_bt_samp   <- sampling_bootstrap(data_resamp_samp,boot_prob,data_signal)
                                      bt_samp_output <- mi_logreg_algorithm(data=data_bt_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                  formula_string=formula_string,model_out=model_out,lr_maxit=lr_maxit,MaxNWts =MaxNWts,
                                                                            pinput=pinput )
                                    bt_samp_output
                                  }
    parallel::stopCluster(cl)
    cat("Boot&Resampling completed..")
  }
  
  # Train-Test
  cl=parallel::makeCluster(testing_cores)
  doParallel::registerDoParallel(cl)
  output_test3=foreach::foreach(j=1:traintest_num,
                                .export=c("sampling_partition","mi_logreg_algorithm","aux_x_log_y"),
                                .packages=c("nnet","caret")) %dopar% {
                                    datatraintestsamp   <- sampling_partition(data,data_signal,partition_trainfrac)
                                    bt_samp_output<- mi_logreg_algorithm(data=datatraintestsamp,signal=signal,response=response,side_variables=side_variables,
                                                                               formula_string=formula_string,model_out=model_out,lr_maxit=lr_maxit,MaxNWts =MaxNWts,
                                                                         pinput=pinput )
                                  bt_samp_output
                                }
  parallel::stopCluster(cl)
  
  print("Over-Fitting-Test completed")
  
  output$bootstrap        <- output_test1
  if (!is.null(side_variables)){ output$resamplingMorph  <- output_test2}
  if (!is.null(side_variables)){ output$bootResampMorph  <- output_test4}
  output$traintest        <- output_test3
  
  output
}
