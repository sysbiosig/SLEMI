#' Probability of pairwise discrimination
#' 
#' INPUT:
#' @param data must be a data.frame object
#' @param signal is a character object that indicates columns of data to be treated as channel's input
#' @param response is a character vector that indicates columns of data to be treated as channel's output
#' @param side_variables is an optional character vector that indicates side variables' columns of data, if NULL no side variables are included
#' @param formula_string is a character object that includes a formula syntax to use in algorithm for capacity calculation
#' @param lr_maxit is the number of iteration of iterative algorithm of logistic regression
#' @param maxNWts is the maximum acceptable number of weights in logistic regression algorithm - the higher, the slower the computations
#' @param model_out is the logical indicating if the calculated logisitc regression model should be included in output list
#' @export
#' @return a list with three elements:
#' \itemize{
#' \item output$prob_matr - confusion matrix of logistic regression predictions
#' \item output$model     - optimal probability distribution
#' }
#' 
prob_discr_pairwise<-function(dataRaw,
                              signal="signal",response="response",side_variables=NULL,
                              formula_string=NULL,
                              output_path=NULL, scale=TRUE,
                              model_out=TRUE,
                              lr_maxit=1000,MaxNWts = 5000){
  
  
  
  print("Procedure starting")
  
  time_start=proc.time()
  
  # checking assumptions
  if (is.null(output_path)) { 
    warning('path is not defined. Graphs and RDS file will not be saved.')
  } else {
    dir.create(output_path,recursive = TRUE)
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
  
  if (!is.null(formula_string)){
    formula=as.formula(formula_string)
  } else {
    formula=as.formula(formula_generator(signal,response, side_variables))
  }
  
  chosen_stim=unique(data[[signal]])
  nstim=length(chosen_stim)
  
  aux_input_checks(data,signal,response,side_variables)
  
  
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
      print("Logistic Regression algorihtm starting")
      lr_model=nnet::multinom(formula,data=dataChosen,na.action=na.omit,maxit=lr_maxit, MaxNWts = MaxNWts)
      model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]$model=lr_model
      
      prob_lr<-data.frame(fitted(lr_model))
      if (length(signal_levels)==2) {prob_lr=cbind(1-prob_lr,prob_lr)}
        
      colnames(prob_lr)<-as.character(signal_levels)
      obs<-dataChosen[[signal]]
      pred<-apply(prob_lr,1,function(x){
          idmax=which.max(x)
          as.character(signal_levels)[idmax]
      })
      model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]$regression<-caret::confusionMatrix(factor(pred,levels=signal_levels),obs)
      
    }
  }
  
  prob_matrix=matrix(0,nstim,nstim)
  for (is in 1:(nstim-1) ){
    for (js in (is+1):nstim){
      w0_acc=c()
      temp_output=model_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]
      w0_acc=sum(temp_output$regression$overall[1])
      prob_matrix[is,js]=w0_acc
    }
  }
  
  prob_matrix=prob_matrix+t(prob_matrix)
  
  for (is in 1:(nstim) ){
      prob_matrix[is,is]=1
  }
  
  row.names(prob_matrix)<-chosen_stim
  colnames(prob_matrix)<-chosen_stim
  
  col2 <- grDevices::colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                             "#FDDBC7",
                             "#67001F", "#B2182B", "#D6604D", "#F4A582",
                             "#FDDBC7","#FFFFFF", "#D1E5F0", "#92C5DE",
                             "#4393C3", "#2166AC", "#053061"))
  
  pdf(paste0(output_path,"/plot_probs_discr.pdf"),height=6,width=6)
    corrplot::corrplot(prob_matrix,type = "upper", method = "pie",cl.lim = c(0.5, 1),col=col2(40), diag=FALSE)
  dev.off()
  
  output=list()
  
  output$prob_matr=prob_matrix
  
  if (model_out){
  output$model=model_output
  }
  
  output
}

