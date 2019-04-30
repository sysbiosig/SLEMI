#' Iterative updating of prior probabilities in logistic regression estimator
#' 
#' Internal, auxillary functions
#'
#' @param prob_lr is a matrix of class probabilties for each observation
#' @param p0 is a numeric vector of prior probabilities used for logistic regression estimation
#' @param cell_id a list of logical vectors indicating class labels of each observation
#' @param signal_levels is a vector of class labels
#' @param cc_maxit is the number of iteration of procedure to be carried out
#' @return A list with components
#' \enumerate{
#'   \item p_opt - a numeric vectors with estimated optimal input probabilities
#'   \item MI_opt -  a numerical value of estimated channel capacity
#' }
#' @keywords internal
func_iterative_logreg_update<-function(prob_lr,p0,cell_id,signal_levels,cc_maxit){
  for (i in 1:cc_maxit){
    C_mc<-sapply(signal_levels,function(x) {
      mc_values=log(prob_lr[[x]][cell_id[[x]] ])
      mean(mc_values[is.finite(mc_values)])
    },simplify=TRUE)
    
    p_opt=exp(C_mc)/sum(exp(C_mc))
    
    #for debugging:
    #MI_opt=sum(C_mc*p_opt-x_log_y(p_opt,p_opt))
    #print(c(p_opt,exp(MI_opt)))  #for debugging
    #print(c(exp(MI_opt)))        #for debugging
    
    prob_ratio=(p0[1]/p0)*(p_opt/p_opt[1])
    temp_val <- prob_lr*t(replicate(nrow(prob_lr),prob_ratio))
    prob_lr <- data.frame(t(apply(temp_val,1,function(x){ x/sum(x) })))
    colnames(prob_lr) <- signal_levels
    p0=p_opt
  }
  MI_opt=sum(C_mc*p_opt-aux_x_log_y(p_opt,p_opt))
  
  out=list(p_opt=p_opt,MI_opt=MI_opt)
  out
}