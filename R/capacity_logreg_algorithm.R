#' Main algorithm to calculate channel capacity by SLEMI approach
#'
#'
#' Additional parameters: lr_maxit and maxNWts are the same as in definition of multinom function from nnet package. An alternative
#' model formula (using formula_string arguments) should be provided if  data are not suitable for description by logistic regression
#' (recommended only for advanced users). It is recommended to conduct estimation by calling capacity_logreg_main.R.
#'
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
#' @param MaxNWts is a maximum acceptable number of weights in logistic regression algorithm. Default is 5000.
#' @param model_out is the logical indicating if the calculated logisitc regression model should be included in output list
#' @export
#' @return a list with three elements:
#' \itemize{
#' \item output$cc         - channel capacity in bits
#' \item output$p_opt      - optimal probability distribution
#' \item output$regression - confusion matrix of logistic regression predictions
#' \item output$model      - nnet object describing logistic regression model (if model_out=TRUE)
#' }
#' @examples 
#' tempdata=data_example1
#' outputCLR1=capacity_logreg_algorithm(data=tempdata, signal="signal", response="response",
#' formula_string = "signal~response")
#'
capacity_logreg_algorithm<-function(data,
                                    signal="signal",response="response",side_variables=NULL,
                                    formula_string=NULL,
                                    model_out=TRUE,
                                    cc_maxit=100,lr_maxit=1000,MaxNWts = 5000){
  
  output<-list()
  
  if (!is.null(formula_string)){
    formula=stats::as.formula(formula_string)
  } else {
    formula=stats::as.formula(func_formula_generator(signal,response, side_variables))
  }
  
  if (is.null(data$train)|is.null(data$test)) {
    
    #cat("\n Main Algorithm starting...")

    # checking assumptions
    func_input_checks(data,signal,response,side_variables)
    
    #transform signal into factor
    data=func_signal_transform(data,signal)
   
    signal_levels<-levels(data[[signal]])
    cell_id=lapply(signal_levels,function(x){
      data[[signal]]==x
    })
    names(cell_id)<-signal_levels
    
    class_num=sapply(cell_id,sum,simplify=TRUE)
    p0=sapply(cell_id,mean,simplify=TRUE)
    if (any(p0==0)) {
      stop('There is no observations of some signals')
    }
    
    
    lr_model=nnet::multinom(formula,data=data,na.action=stats::na.omit,maxit=lr_maxit, MaxNWts = MaxNWts,trace=FALSE)
    
    prob_lr<-data.frame(stats::fitted(lr_model))
    if (length(signal_levels)==2) {prob_lr=cbind(1-prob_lr,prob_lr)}
    
    colnames(prob_lr)<-as.character(signal_levels)
    obs<-data[[signal]]
    pred<-apply(prob_lr,1,function(x){
      idmax=which.max(x)
      as.character(signal_levels)[idmax]
    })
    output$regression<-caret::confusionMatrix(factor(pred,levels=signal_levels),obs)
    rm(obs,pred)
    
    #Debugging:
    #print("Iterative algorihtm starting")
    
    tmp_iterative_output=func_iterative_logreg_update(prob_lr,p0,cell_id,signal_levels,cc_maxit)
    
    #Debugging:
    #print("Iterative algorihtm complete")
    
    if (model_out) {
      output$model = lr_model
      #temp_z <- coef(output$model)/(summary(output$model)$standard.errors)
      #output$model_pv <- (1 - pnorm(abs(temp_z), 0, 1)) * 2
    }
    output$p_opt  <- tmp_iterative_output$p_opt
    output$cc     <- log2(exp(tmp_iterative_output$MI_opt))
    
    #  cat("... Main algorithm completed.")

  } else {
    
    data_train<-data$train
    data_test <-data$test
    
    # checking assumptions
    func_input_checks(data_train,signal,response,side_variables)
    func_input_checks(data_test,signal,response,side_variables)
    
    #transform signal into factor
    data_train=func_signal_transform(data_train,signal)
    data_test=func_signal_transform(data_test,signal)
    
    signal_levels<-levels(data_train[[signal]])
    
    cell_id_train=lapply(signal_levels,function(x){
      data_train[[signal]]==x
    })
    names(cell_id_train)<-signal_levels
    cell_id_test=lapply(signal_levels,function(x){
      data_test[[signal]]==x
    })
    names(cell_id_test)<-signal_levels
    
    class_num=sapply(cell_id_train,sum,simplify=TRUE)
    p0=sapply(cell_id_train,mean,simplify=TRUE)
    if (any(p0==0)) {
      stop('There is no observations of some signals')
    }
    
    lr_model=nnet::multinom(formula,data=data_train,na.action=stats::na.omit,maxit=lr_maxit, MaxNWts = MaxNWts,trace=FALSE)
    
    prob_lr_train<-data.frame(stats::fitted(lr_model))
    if (length(signal_levels)==2) {prob_lr_train=cbind(1-prob_lr_train,prob_lr_train)}
    prob_lr_test<-data.frame(stats::predict(lr_model,newdata=data_test,type="prob"))
    if (length(signal_levels)==2) {prob_lr_test=cbind(1-prob_lr_test,prob_lr_test)}
    
    colnames(prob_lr_train)<-as.character(signal_levels)
    colnames(prob_lr_test)<-as.character(signal_levels)
    
    obs_train<-data_train[[signal]]
    obs_test<-data_test[[signal]]
    
    pred_train<-apply(prob_lr_train,1,function(x){
      idmax=which.max(x)
      as.character(signal_levels)[idmax]
    })
    pred_test<-apply(prob_lr_test,1,function(x){
      idmax=which.max(x)
      as.character(signal_levels)[idmax]
    })
    
    output$regression<-list()
    output$regression$test<-caret::confusionMatrix(factor(pred_test,levels=signal_levels),obs_test)
    output$regression$train<-caret::confusionMatrix(factor(pred_train,levels=signal_levels),obs_train)
    
    tmp_iterative_output=func_iterative_logreg_update(prob_lr_test,p0,cell_id_test,signal_levels,cc_maxit)
  
    output$p_opt  <- tmp_iterative_output$p_opt
    output$cc     <- log2(exp(tmp_iterative_output$MI_opt))
    
    if (model_out){
      #output$data   <- data # for debugging
      output$model <- lr_model
      #temp_z <- coef(output$model)/(summary(output$model)$standard.errors)
      #output$model_pv <- (1 - pnorm(abs(temp_z), 0, 1)) * 2
    }
  }

  
  output
}

