#' Plotting output of capacity estimation and additional exploratory graphs.
#' 
#' INPUT:
#' @param data is a data.frame object
#' @param signal is a character object that indicates columns of data that should be treated as channel's input
#' @param response is a character vector that indicates columns of data that should be treated as channel's output
#' @param side_variables is a character vector that indicates side variables' columns of data
#' @param cc_output a list that is a standard output of capacity_logreg_algorithm function
#' @param path character giving the directory, where graphs should be saved
#' @param height integer indicating the height of a single plot
#' @param width integer indicating the width of a single plot
#' @return A list with ggplot or gtable object. Each plot is also saved in `output_path` directory in separate pdfs files which include:
#' \itemize{
#' \item MainPlot.pdf - a simple summary plot with basic distribution visualization and capacity estimate
#' \item MainPlot_full.pdf - a comprehensive summary plot with distribution visualization and capacity estimate
#' \item capacity.pdf - a diagram presenting the capacity estimates
#' \item io_relation.pdf - a graph with input-output relation
#' \item kdensities.pdf - kernel density estimator of data distribution
#' \item histograms.pdf - histograms of data
#' \item boxplots.pdf - boxplots of data
#' \item violin.pdf - violin plots of data
#' }
#' @keywords internal
#' 
output_graphs_main<-function(data,signal,response,side_variables,cc_output,
                        output_path,height=4,width=6){
  
  data=func_signal_transform(data,signal)
  warn_ind=0

  response_factor=factor(response)
  response_class=levels(response_factor)
  if (length(response_class)>3){
    response_class=response_class[seq(from=1,to=length(response_class),length.out = 3)]
    response=response[response%in%response_class]
    message("There are more than 3 dimensions in output. Restricting plots to 3 (minimal, middle and maximal).")
  }

  plot_box=try(capacity_output_graph_boxplots(data,signal,response,output_path,height=height,width=width),silent=FALSE)
  plot_violin=try(capacity_output_graph_violinMean(data,signal,response,output_path,height=height,width=width),silent=FALSE)
  plot_boxSideVar=try(capacity_output_graph_boxplotsSideVar(data,signal,side_variables,output_path,height=height,width=width) ,silent=FALSE)
  plot_capacity=try(capacity_output_graph_capacity(cc_output,output_path,height=height,width=width),silent=FALSE)

  if (any(class(plot_box)=="try-error")){
    plot_box=grid::textGrob("")
    warn_ind=1
  }
  if (any(class(plot_violin)=="try-error")){
    plot_violin=grid::textGrob("")
    warn_ind=1
  } 
  if (any(class(plot_boxSideVar)=="try-error")){
    plot_boxSideVar=grid::textGrob("")
    warn_ind=1
  }
  if (any(class(plot_capacity)=="try-error")){
    plot_capacity=grid::textGrob("")
    warn_ind=1
  }

  if (is.null(cc_output$testing)){
    plot_main_simp=try(gridExtra::grid.arrange(plot_violin,plot_box,plot_capacity,
                                          layout_matrix=rbind(c(1,2),
                                                              c(1,2),
                                                              c(1,2),
                                                              c(3,3),
                                                              c(3,3)
                                          )),silent=FALSE)
    if(!is.null(output_path)){
      try(ggplot2::ggsave(plot_main_simp,file=paste(output_path,'MainPlot.pdf',sep=""),
        height=2*height,width=2*width,limitsize = FALSE),silent=FALSE)
    }

    if (any(class(plot_main_simp)=="try-error")){
      warn_ind=1
    }

  } else {
    plot_main_simp=try(gridExtra::grid.arrange(plot_violin,plot_box,plot_capacity,
                                          layout_matrix=rbind(c(1,2),
                                                              c(3,3),
                                                              c(3,3),
                                                              c(3,3)
                                          )),silent=FALSE)
    if(!is.null(output_path)){
      try(ggplot2::ggsave(plot_main_simp,file=paste(output_path,'MainPlot.pdf',sep=""),
        height=3*height,width=2*width,limitsize = FALSE),silent=FALSE)
    }

    if (any(class(plot_main_simp)=="try-error")){
      warn_ind=1
    }
  }
  
  if (warn_ind==1){
    warning("At least one diagnostic plots has not been created. Check input data.")
  }
  
  
  graphOutput=list(plot_box,plot_boxSideVar,plot_capacity,plot_main_simp)
}