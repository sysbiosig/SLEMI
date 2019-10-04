#' Testing procedures for estimation of channel capacity
#' 
#' Diagnostic procedures that allows to compute the uncertainity of estimation of channel capacity by SLEMI approach. Two main procedures are implemented:
#' bootstrap, which execute estimation with using a fraction of data and overfitting test, which dividies data into two parts: training and testing. Each of them
#' is repeated specified number of times to obtain a distribution of our estimators. It is recommended to conduct estimation by calling capacity_logreg_main.R.
#'
#' If side variables are added within the analysis (side_variables is not NULL), two additional procedures are carried out:
#' reshuffling test and reshuffling with boostrap test, which are based on permutation of side variables values within the dataset.  
#' Additional parameters: lr_maxit and maxNWts are the same as in definition of multinom function from nnet package. An alternative
#' model formula (using formula_string arguments) should be provided if  data are not suitable for description by logistic regression
#' (recommended only for advanced users).
#'
#' @section References:
#' [1] Jetka T, Nienaltowski K, Winarski T, Blonski S, Komorowski M,  
#' Information-theoretic analysis of multivariate single-cell signaling responses using SLEMI,
#' \emph{PLoS Comput Biol}, 15(7): e1007132, 2019, https://doi.org/10.1371/journal.pcbi.1007132.
#'
#' @param data must be a data.frame object. Cannot contain NA values.
#' @param signal is a character object with names of columns of dataRaw to be treated as channel's input.
#' @param response is a character vector with names of columns of dataRaw  to be treated as channel's output
#' @param side_variables (optional) is a character vector that indicates side variables' columns of data, if NULL no side variables are included
#' @param formula_string (optional) is a character object that includes a formula syntax to use in logistic regression model. 
#' If NULL, a standard additive model of response variables is assumed. Only for advanced users.
#' @param cc_maxit is the number of iteration of iterative optimisation of the algorithm to esimate channel capacity. Default is 100.
#' @param lr_maxit is a maximum number of iteration of fitting algorithm of logistic regression. Default is 1000.
#' @param maxNWts is a maximum acceptable number of weights in logistic regression algorithm. Default is 5000.
#' @param TestingSeed is the seed for random number generator used in testing procedures
#' @param boot_num is the number of bootstrap tests to be performed. Default is 10, but it is recommended to use at least 50 for reliable estimates.
#' @param boot_prob is the proportion of initial size of data to be used in bootstrap. Default is 0.8.
#' @param testing_cores - number of cores to be used in parallel computing (via doParallel package)
#' @param sidevar_num is the number of re-shuffling tests of side variables to be performed. Default is 10, but it is recommended to use at least 50 for reliable estimates.
#' @param traintest_num is the number of overfitting tests to be performed. Default is 10, but it is recommended to use at least 50 for reliable estimates.
#' @param partition_trainfrac is the fraction of data to be used as a training dataset. Default is 0.6.
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
#' ## Please set boot_num and traintest_num with larger numbers 
#' ## for a more reliable testing
#' tempdata=data_example1
#' outputCLR1_testing=capacity_logreg_testing(data=tempdata,
#' signal="signal", response="response",cc_maxit=10,
#' TestingSeed=11111, boot_num=1,boot_prob=0.8,testing_cores=1,
#' traintest_num=1,partition_trainfrac=0.6)
#'
capacity_logreg_testing<-function(data,signal="signal",response="response",side_variables=NULL,
                                          cc_maxit=100,lr_maxit=1000,MaxNWts = 5000,
                                          formula_string=NULL,
                                          TestingSeed=1234,testing_cores=1,
                                          boot_num=10,boot_prob=0.8,
                                          sidevar_num=10,
                                          traintest_num=10,partition_trainfrac=0.6){
  
  output=list()
  data_signal=data[[signal]]
  
  set.seed(TestingSeed)
  message("Testing procedures starting with ",testing_cores," core(s)..")
  
    if (testing_cores==1){  
        `%do%`<-foreach::`%do%`
        output_test1<-foreach::foreach(j=1:boot_num,
                                   .export=c("sampling_bootstrap","capacity_logreg_algorithm","aux_x_log_y","func_formula_generator",
                                    "func_input_checks","func_signal_transform","func_iterative_logreg_update"),
                                   .packages=c("nnet","caret")) %do% {
                                        data_bt_samp   <- sampling_bootstrap(data,boot_prob,data_signal)
                                        bt_samp_output <- capacity_logreg_algorithm(data=data_bt_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                              formula_string=formula_string,model_out=FALSE,cc_maxit=cc_maxit,
                                                                                              lr_maxit=lr_maxit,MaxNWts =MaxNWts )
                                        bt_samp_output
                                    }
      
        if (!is.null(side_variables)){
            output_test2=foreach::foreach(j=1:sidevar_num,
                                            .export=c("sampling_shuffle","capacity_logreg_algorithm","aux_x_log_y","func_formula_generator",
                                            "func_input_checks","func_signal_transform","func_iterative_logreg_update"),
                                            .packages=c("nnet","caret")) %do% {
                                                data_resamp_samp   <- sampling_shuffle(data,side_variables)
                                                bt_samp_output <- capacity_logreg_algorithm(data=data_resamp_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                            formula_string=formula_string,model_out=FALSE,cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts )
                                                bt_samp_output
                                            }
            output_test4=foreach::foreach(j=1:sidevar_num,
                                        .export=c("sampling_bootstrap","sampling_shuffle","capacity_logreg_algorithm","aux_x_log_y","func_formula_generator",
                                            "func_input_checks","func_signal_transform","func_iterative_logreg_update"),
                                        .packages=c("nnet","caret")) %do% {
                                            data_resamp_samp   <- sampling_shuffle(data,side_variables)
                                            data_bt_samp   <- sampling_bootstrap(data_resamp_samp,boot_prob,data_signal)
                                            bt_samp_output <- capacity_logreg_algorithm(data=data_bt_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                        formula_string=formula_string,model_out=FALSE,cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts )
                                            bt_samp_output
                                        }
        }
  
        output_test3=foreach::foreach(j=1:traintest_num,
                                .export=c("sampling_partition","capacity_logreg_algorithm","aux_x_log_y","func_formula_generator",
                                  "func_input_checks","func_signal_transform","func_iterative_logreg_update"),
                                .packages=c("nnet","caret")) %do% {
                                        datatraintestsamp   <- sampling_partition(data,data_signal,partition_trainfrac)
                                        bt_samp_output<- capacity_logreg_algorithm(data=datatraintestsamp,signal=signal,response=response,side_variables=side_variables,
                                                                                 formula_string=formula_string,model_out=FALSE,cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts )
                                        bt_samp_output
                                    }

        output$bootstrap        <- output_test1
        if (!is.null(side_variables)){ output$reshuffling_sideVar  <- output_test2}
        if (!is.null(side_variables)){ output$bootstrap_Reshuffling_sideVar  <- output_test4}
        output$traintest        <- output_test3
    } else if (testing_cores>1) {
        `%dopar%`<-foreach::`%dopar%`
        
        #Bootstrap
        #message(" Bootstrap starting..")
        cl=parallel::makeCluster(testing_cores)
        doParallel::registerDoParallel(cl)
        output_test1<-foreach::foreach(j=1:boot_num,
                                       .export=c("sampling_bootstrap","capacity_logreg_algorithm","aux_x_log_y","func_formula_generator",
                                        "func_input_checks","func_signal_transform","func_iterative_logreg_update"),
                                       .packages=c("nnet","caret")) %dopar% {
                                            data_bt_samp   <- sampling_bootstrap(data,boot_prob,data_signal)
                                            bt_samp_output <- capacity_logreg_algorithm(data=data_bt_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                                formula_string=formula_string,model_out=FALSE,cc_maxit=cc_maxit,
                                                                                                lr_maxit=lr_maxit,MaxNWts =MaxNWts )
                                          bt_samp_output
                                        }
        #parallel::stopCluster(cl)
        #message("... completed")
        
        if (!is.null(side_variables)){
        #Resampling
        #message(" Reshuffling starting..")
        #cl=parallel::makeCluster(testing_cores)
        #doParallel::registerDoParallel(cl)
            output_test2=foreach::foreach(j=1:sidevar_num,
                                        .export=c("sampling_shuffle","capacity_logreg_algorithm","aux_x_log_y","func_formula_generator",
                                          "func_input_checks","func_signal_transform","func_iterative_logreg_update"),
                                        .packages=c("nnet","caret")) %dopar% {
                                            data_resamp_samp   <- sampling_shuffle(data,side_variables)
                                            bt_samp_output <- capacity_logreg_algorithm(data=data_resamp_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                      formula_string=formula_string,model_out=FALSE,cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts )
                                            bt_samp_output
                                        }
        #parallel::stopCluster(cl)
        #message("... completed 1 ...")
        
        #Bootstrap&Resampling
        #cl=parallel::makeCluster(testing_cores)
        #doParallel::registerDoParallel(cl)
            output_test4=foreach::foreach(j=1:sidevar_num,
                                      .export=c("sampling_bootstrap","sampling_shuffle","capacity_logreg_algorithm","aux_x_log_y","func_formula_generator",
                                        "func_input_checks","func_signal_transform","func_iterative_logreg_update"),
                                      .packages=c("nnet","caret")) %dopar% {
                                            data_resamp_samp   <- sampling_shuffle(data,side_variables)
                                            data_bt_samp   <- sampling_bootstrap(data_resamp_samp,boot_prob,data_signal)
                                            bt_samp_output <- capacity_logreg_algorithm(data=data_bt_samp,signal=signal,response=response,side_variables=side_variables,
                                                                                        formula_string=formula_string,model_out=FALSE,cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts )
                                            bt_samp_output
                                      }
        #parallel::stopCluster(cl)
        #message("completed 2")
        }
        
        # Train-Test
        #message(" Over-fitting starting..")
        # cl=parallel::makeCluster(testing_cores)
        # doParallel::registerDoParallel(cl)
        output_test3=foreach::foreach(j=1:traintest_num,
                                      .export=c("sampling_partition","capacity_logreg_algorithm","aux_x_log_y","func_formula_generator",
                                        "func_input_checks","func_signal_transform","func_iterative_logreg_update"),
                                      .packages=c("nnet","caret")) %dopar% {
                                            datatraintestsamp   <- sampling_partition(data,data_signal,partition_trainfrac)
                                            bt_samp_output<- capacity_logreg_algorithm(data=datatraintestsamp,signal=signal,response=response,side_variables=side_variables,
                                                                                       formula_string=formula_string,model_out=FALSE,cc_maxit=cc_maxit,lr_maxit=lr_maxit,MaxNWts =MaxNWts )
                                            bt_samp_output
                                        }
        parallel::stopCluster(cl)
        #message("... completed")
        
        output$bootstrap        <- output_test1
        if (!is.null(side_variables)){ output$reshuffling_sideVar  <- output_test2}
        if (!is.null(side_variables)){ output$bootstrap_Reshuffling_sideVar  <- output_test4}
        output$traintest        <- output_test3
    }
      
    output
}
