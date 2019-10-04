#' Calculates Probability of pairwise discrimination
#' 
#' Estimates probabilities of correct discrimination (PCDs) between each pair of input/signal values using a logistic regression model.
#'
#' 
#' In order to estimate PCDs, for a given pair of input values \eqn{x_i} and \eqn{x_j}, we propose to fit a logistic regression model
#' using response data corresponding to the two considered inputs, i.e. \eqn{y^l_u}, for \eqn{l\in\{i,j\}} and \eqn{u} ranging from 
#' 1 to \eqn{n_l}. 
#' To ensure that both inputs have equal contribution to the calculated discriminability, equal probabilities should be assigned,
#' \eqn{P(X) = (P(x_i),P(x_j))=(1/2,1/2)}. Once the regression model is fitted, probability of assigning a given cellular response, 
#' \eqn{y},
#' to the correct input value is estimated as
#' \deqn{\max \{ \hat{P}_{lr}(x_i|Y=y;P(X)), \hat{P}_{lr}(x_j|Y=y;P(X))\}.}
#' Note that \eqn{P(x_j|Y=y)=1-P(x_i|Y=y)} as well as \eqn{\hat{P}_{lr}(x_j|Y=y;P(X))=1-\hat{P}_{lr}(x_i|Y=y;P(X))}
#' The average of the above probabilities over all observations \eqn{y^i_l} yields PCDs
#' \deqn{PCD_{x_i,x_j}=\frac{1}{2}\frac{1}{n_i}\sum_{l=1}^{n_i}\max\{ \hat{P}_{lr}(x_i|Y=y_i^l;P(X)),\hat{P}_{lr}(x_i^l|Y=y;P(X))\} + }
#' \deqn{ \frac{1}{2}  \frac{1}{n_j} \sum_{l=1}^{n_j} \max \{ \hat{P}_{lr}(x_i|Y=y_j^l;P(X)), \hat{P}_{lr}(x_j|Y=y_j^l;P(X))\}.}
#' 
#' Additional parameters: lr_maxit and maxNWts are the same as in definition of multinom function from nnet package. An alternative
#' model formula (using formula_string arguments) should be provided if  data are not suitable for description by logistic regression
#' (recommended only for advanced users). Preliminary scaling of  data (argument scale) should be used similarly as in other 
#' data-driven approaches, e.g. if response variables are comparable, scaling (scale=FALSE) can be omitted, while if they represent 
#' different phenomenon (varying by units and/or magnitude) scaling is recommended.
#'
#' @section References:
#' [1] Jetka T, Nienaltowski K, Winarski T, Blonski S, Komorowski M,  
#' Information-theoretic analysis of multivariate single-cell signaling responses using SLEMI,
#' \emph{PLoS Comput Biol}, 15(7): e1007132, 2019, https://doi.org/10.1371/journal.pcbi.1007132.
#'
#' @param dataRaw must be a data.frame object
#' @param signal is a character object with names of columns of dataRaw to be treated as channel's input.
#' @param response is a character vector with names of columns of dataRaw  to be treated as channel's output
#' @param side_variables (optional) is a character vector that indicates side variables' columns of data, if NULL no side variables are included
#' @param formula_string (optional) is a character object that includes a formula syntax to use in logistic regression model. 
#' If NULL, a standard additive model of response variables is assumed. Only for advanced users.
#' @param output_path is a directory where a pie chart with calculated probabilities will be saved. If NULL, the graph will not be created.
#' @param scale is a logical indicating if the response variables should be scaled and centered before fitting logistic regression
#' @param lr_maxit is a maximum number of iteration of fitting algorithm of logistic regression. Default is 1000.
#' @param MaxNWts is a maximum acceptable number of weights in logistic regression algorithm. Default is 5000.
#' @param diagnostics is a logical indicating if details of logistic regression fitting should be included in output list
#' @export
#' @return a list with two elements:
#' \itemize{
#' \item output$prob_matr - a \eqn{n\times n} matrix, where \eqn{n} is the number of inputs, with probabilities of correct 
#' discirimination between pairs of input values. 
#' \item output$diagnostics     - (if diagnostics=TRUE) a list correspondning to logistic regression models fitted for each 
#' pair of input values. Each element consists of three sub-elements: 1) nnet_model - nnet object summarising logistic regression model; 
#' 2) prob_lr - probabilities of assignment obtained from logistic regression model; 
#' 3) confusion_matrix - confusion matrix of classificator.
#' }
#' @examples 
#' ## Calculate probabilities of discrimination for toy dataset
#' temp_data=data_example1
#' output=prob_discr_pairwise(dataRaw=data_example1,
#'                    signal = "signal",
#'                    response = "response")
#' ## Calculate probabilities of discrimination for nfkb dataset
#'  it=21 # choose from 0, 3, 6, ..., 120 for measurements at other time points
#'  output=prob_discr_pairwise(dataRaw=data_nfkb,
#'                             signal = "signal",
#'                            response = paste0("response_",it))
#'
prob_discr_pairwise<-function(dataRaw,
                              signal="input",response=NULL,side_variables=NULL,
                              formula_string=NULL,
                              output_path=NULL, scale=TRUE,
                              lr_maxit=1000,MaxNWts = 5000,diagnostics=TRUE){
  
  
  
  message("Estimating pairwise probabilities of discrimination...")
  
  if (is.null(response)){
    response=paste0("output_",1:(ncol(dataRaw)-1) )  
  }
  
  time_start=proc.time()
  dataRaw=as.data.frame(dataRaw)

  # checking assumptions
  if (is.null(output_path)) { 
    message('Path is not defined. Graphs and RDS file will not be saved.')
  } else {
    dir.create(output_path,recursive = TRUE)
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
    message(" There are NA in observations - removing..")
    data0=data0[!apply(data0,1,function(x) any(is.na(x)) ),]
  }
  

  data0=func_signal_transform(data0,signal)
  tempcolnames=colnames(data0)
  tempsignal=data.frame(data0[,(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )])
  colnames(tempsignal)<-tempcolnames[(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )]
  data0=data.frame(data0[,!(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )])
  colnames(data0)<-tempcolnames[!(tempcolnames%in%c(signal,paste(signal,"_RAW",sep="") ) )]
  
 # message("  Preprocessing started...")
  
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
  #cat(" completed")
  
  if (!is.null(formula_string)){
    formula=stats::as.formula(formula_string)
  } else {
    formula=stats::as.formula(func_formula_generator(signal,response, side_variables))
  }
  
  chosen_stim=unique(data[[signal]])
  nstim=length(chosen_stim)
  pinput=c(0.5,0.5)

  func_input_checks(data,signal,response,side_variables)
  
  temp_num_pairs=(nstim^2-nstim)*0.5
  message("Fitting logistic regression models for ",temp_num_pairs," pairs")
  model_output=list()
  # 11 Estimate classificator
  for (is in 1:(nstim-1) ){
    for (js in (is+1):nstim){

      capacity_chosen_stim=chosen_stim[c(is,js)]
      model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]=list()
          
      dataChosen=data[data[[signal]]%in% capacity_chosen_stim,]
      dataChosen[[signal]]=factor(as.character(dataChosen[[signal]]))
      
      signal_levels<-levels(dataChosen[[signal]])
      cell_id=lapply(signal_levels,function(x){
        dataChosen[[signal]]==x
      })
      names(cell_id)<-signal_levels
        
      class_num=sapply(cell_id,sum,simplify=TRUE)
      p0=sapply(cell_id,mean,simplify=TRUE)
      if (any(p0==0)) {
        stop('There is no observations of some signals')
      }

      #Debugging:
      lr_model=nnet::multinom(formula,data=dataChosen,na.action=stats::na.omit,maxit=lr_maxit, MaxNWts = MaxNWts,trace=FALSE)
      model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]$nnet_model=lr_model
      
      prob_lr<-data.frame(stats::fitted(lr_model))
      if (length(signal_levels)==2) {prob_lr=cbind(1-prob_lr,prob_lr)}
        
      prob_ratio=(p0[1]/p0)*(pinput/pinput[1])
      temp_val <- prob_lr*t(replicate(nrow(prob_lr),prob_ratio))
      prob_lr <- data.frame(t(apply(temp_val,1,function(x){ x/sum(x) })))
      colnames(prob_lr) <- signal_levels

      obs<-dataChosen[[signal]]
      pred<-apply(prob_lr,1,function(x){
          idmax=which.max(x)
          as.character(signal_levels)[idmax]
      })

      prob_lr=data.frame(signal=dataChosen[[signal]],prob_lr)
      model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]$prob_lr<-prob_lr
      model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]$confusion_matrix<-caret::confusionMatrix(factor(pred,levels=signal_levels),obs)
      
    }
  }
  #cat("completed...")

  prob_matrix=matrix(0,nstim,nstim)
  for (is in 1:(nstim-1) ){
    for (js in (is+1):nstim){
      w0_acc=c()
      #temp_output=model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]
      #w0_acc=sum(temp_output$regression$overall[1])
      #prob_matrix[is,js]=w0_acc
      temp_probs=model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]$prob_lr
      prob_matrix[is,js]=mean(c(by(temp_probs,temp_probs$signal,function(x){
        mean(apply(x[,-1],1,max))
        })))
    }
  }
  
  prob_matrix=prob_matrix+t(prob_matrix)
  
  for (is in 1:(nstim) ){
      prob_matrix[is,is]=1
  }
  
  row.names(prob_matrix)<-chosen_stim
  colnames(prob_matrix)<-chosen_stim
  
  col2 <- grDevices::colorRampPalette(c(rep("#FFFFFF",16), "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))

  #prob_matrix[2,3]=0
  #prob_matrix=prob_matrix*0.5
  #prob_matrix[prob_matrix<0.5]=0.5
  #print(prob_matrix)

min_val=floor(min(prob_matrix)*10)/10
if (min_val<0.5){min_val=-1}


if (!is.null(output_path)){                       
  grDevices::pdf(paste0(output_path,"/plot_probs_discr.pdf"),height=0.75*length(chosen_stim),width=0.75*length(chosen_stim))
    corrplot::corrplot(prob_matrix,type = "upper", method = "pie",
      cl.lim = c( min(c(min_val,0.5)) , 1),col=col2(100), diag=FALSE)
  grDevices::dev.off()
   message("Graph with probabilities created in",output_path) 
  }
                     
                       
  for (is in 1:(nstim) ){
      prob_matrix[is,is]=NA
  }
  
  output=list()
  
  output$prob_matr=prob_matrix
  
  if (diagnostics){
    output$diagnostics=model_output
  }
  
  output
}

