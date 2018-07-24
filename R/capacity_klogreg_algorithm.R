#' Algorithm for channel capacity calculation using glmnet algorithm
#' 
#' @param dataMatrix is a numeric matrix with columns treated as explanatory variables in 
#' logistic regression algorithm (only if glmnet_algorithm=TRUE)
#' @param dataSignal is a vector of factors indicationg channel's input
#' @param cv_core_num is the number of cores to use in parallel computing of glmnet package 
#' @param lambda_num is the lambda parameter of glmnet package
#' @param cc_maxit is the number of iteration of iterative algorithm of finding maximum CC
#' @param model_out is the logical indicating if the calculated logisitc regression model should be included in output list
#'
#' @return a list with three elements:
#' \itemize{
#' \item output$regression - confusion matrix of logistic regression predictions
#' \item output$cc         - channel capacity in bits
#' \item output$p_opt      - optimal probability distribution
#' \item output$model_out      - glmnet object describing logistic regression model (if model_out=TRUE)
#' }
capacity_klogreg_algorithm<-function(dataMatrix,
                                     dataSignal,
                                                  model_out=TRUE,
                                                  cc_maxit=100,
                                                  lambda_num=10,
                                                  cv_core_num=1){
  
  output<-list()
  
  data=dataMatrix
  
  if (!(is.list(data)|is.list(data)) ) {
    
    
    # checking assumptions
    if (!is.matrix(data)) {
      stop('data is not in matrix format')
    }
    if ( length(dataSignal)==0 ) {
      stop('There is no column described as signal in data')
    }
    
    
    signal_levels<-levels(dataSignal)
    cell_id=lapply(signal_levels,function(x){
      dataSignal==x
    })
    names(cell_id)<-signal_levels
    
    class_num=sapply(cell_id,sum,simplify=TRUE)
    p0=sapply(cell_id,mean,simplify=TRUE)
    if (any(p0==0)) {
      stop('There is no observations of some signals')
    }
    
    if (!is.null(cv_core_num)){
      cl<-parallel::makeCluster(cv_core_num)
      doParallel::registerDoParallel(cl)
      if (length(signal_levels)==2) {lr_cvfit = glmnet::cv.glmnet(dataMatrix,  dataSignal, family = "binomial", type.measure = "mae", alpha = 0.5, nlambda = lambda_num, parallel=TRUE)}
      if (length(signal_levels)>2)  {lr_cvfit = glmnet::cv.glmnet(dataMatrix, dataSignal, family = "multinomial", type.multinomial = "grouped", type.measure = "mae", alpha = 0.5, nlambda = lambda_num, parallel=TRUE)}  
      parallel::stopCluster(cl)
    } else {
      if (length(signal_levels)==2) {lr_cvfit = glmnet::cv.glmnet(dataMatrix,  dataSignal, family = "binomial", type.measure = "mae", alpha = 0.5, nlambda = lambda_num, parallel=FALSE)}
      if (length(signal_levels)>2)  {lr_cvfit = glmnet::cv.glmnet(dataMatrix, dataSignal, family = "multinomial", type.multinomial = "grouped", type.measure = "mae", alpha = 0.5, nlambda = lambda_num, parallel=FALSE)}  
    }
    
    prob_lr<-data.frame(predict(lr_cvfit, newx=dataMatrix, s = "lambda.1se" ,type="response"))
    
    if (length(signal_levels)==2) {prob_lr=cbind(1-prob_lr,prob_lr)}
    
    colnames(prob_lr)<-as.character(signal_levels)
    obs<-dataSignal
    pred<-apply(prob_lr,1,function(x){
      idmax=which.max(x)
      as.character(signal_levels)[idmax]
    })
    output$regression<-caret::confusionMatrix(factor(pred,levels=signal_levels),obs)
    rm(obs,pred)
    
    #Debugging:
    print("Iterative algorihtm starting")
    
    tmp_iterative_output=aux_iterative_logreg_update(prob_lr,p0,cell_id,signal_levels,cc_maxit)
    
    #Debugging:
    print("Iterative algorihtm complete")
    
    if (model_out) {
      output$model = lr_cvfit
    }
    output$p_opt  <- tmp_iterative_output$p_opt
    output$cc     <- log2(exp(tmp_iterative_output$MI_opt))
    
  } else {
    
    data_train<-data$train
    data_test <-data$test
    signal_train<-dataSignal$train
    signal_test<-dataSignal$test
    
    # checking assumptions
    if (!is.matrix(data_test)|!is.matrix(data_train)) {
      stop('data is not in matrix format')
    }
    if ( length(signal_train)==0|length(signal_test)==0 ) {
      stop('There is no column described as signal in data')
    }
    
    signal_levels<-levels(signal_train)
    
    cell_id_train=lapply(signal_levels,function(x){
      signal_train==x
    })
    names(cell_id_train)<-signal_levels
    cell_id_test=lapply(signal_levels,function(x){
      signal_test==x
    })
    names(cell_id_test)<-signal_levels
    
    class_num=sapply(cell_id_train,sum,simplify=TRUE)
    p0=sapply(cell_id_train,mean,simplify=TRUE)
    if (any(p0==0)) {
      stop('There is no observations of some signals')
    }
    
    if (!is.null(cv_core_num)){
      cl<-parallel::makeCluster(cv_core_num)
      doParallel::registerDoParallel(cl)
      if (length(signal_levels)==2) {lr_cvfit = glmnet::cv.glmnet(data_train,   signal_train, family = "binomial", type.measure = "mae", alpha = 0.5, nlambda = lambda_num, parallel=TRUE)}
      if (length(signal_levels)>2)  {lr_cvfit = glmnet::cv.glmnet(data_train,  signal_train, family = "multinomial", type.multinomial = "grouped", type.measure = "mae", alpha = 0.5, nlambda = lambda_num, parallel=TRUE)}  
      parallel::stopCluster(cl)
    } else {
      if (length(signal_levels)==2) {lr_cvfit = glmnet::cv.glmnet(data_train,   signal_train, family = "binomial", type.measure = "mae", alpha = 0.5, nlambda = lambda_num, parallel=FALSE)}
      if (length(signal_levels)>2)  {lr_cvfit = glmnet::cv.glmnet(data_train, signal_train, family = "multinomial", type.multinomial = "grouped", type.measure = "mae", alpha = 0.5, nlambda = lambda_num, parallel=FALSE)}  
    }
    
    prob_lr_train<-data.frame(predict(lr_cvfit, newx=data_train, s = "lambda.1se" ,type="response"))
    if (length(signal_levels)==2) {prob_lr_train=cbind(1-prob_lr_train,prob_lr_train)}
    prob_lr_test<-data.frame(predict(lr_cvfit, newx=data_test, s = "lambda.1se" ,type="response"))
    if (length(signal_levels)==2) {prob_lr_test=cbind(1-prob_lr_test,prob_lr_test)}
    
    colnames(prob_lr_train)<-as.character(signal_levels)
    colnames(prob_lr_test)<-as.character(signal_levels)
    
    obs_train<-signal_train
    obs_test<-signal_test
    
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
  
    
    tmp_iterative_output=aux_iterative_logreg_update(prob_lr_test,p0,cell_id_test,signal_levels,cc_maxit)
    
    output$p_opt  <- tmp_iterative_output$p_opt
    output$cc     <- log2(exp(tmp_iterative_output$MI_opt))
    
    if (model_out){
     # output$data   <- data # for debugging
      output$model <- lr_cvfit
    }
  }
  
  #Debugging:
  print("Main algorihtm complete")
  
  output
}