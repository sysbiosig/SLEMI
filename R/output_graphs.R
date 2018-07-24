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
#' @keywords data.frame
#' @examples 
#' 
#' 
capacity_output_graphs<-function(data,signal,response,side_variables,cc_output,
                        output_path,height=4,width=6){
  
  data=aux_signal_transform(data,signal)
  
  plot_io=try(capacity_output_graph_io(data,signal,response,output_path,height=height,width=width),silent=TRUE)
  plot_box=try(capacity_output_graph_boxplots(data,signal,response,output_path,height=height,width=width),silent=FALSE)
  plot_violin=try(capacity_output_graph_violinMean(data,signal,response,output_path,height=height,width=width),silent=FALSE)
  plot_hist=try(capacity_output_graph_histograms(data,signal,response,output_path,height=height,width=width),silent=FALSE)
  plot_boxSideVar=try(capacity_output_graph_boxplotsSideVar(data,signal,side_variables,output_path,height=height,width=width) ,silent=FALSE)
  plot_capacity=try(capacity_output_graph_capacity(cc_output,output_path,height=height,width=width),silent=FALSE)
  plot_density=try(capacity_output_graph_densities(data,signal,response,output_path,height=height,width=width),silent=FALSE)
  
  if (!is.null(side_variables)){
    if (is.null(cc_output$testing)){
      plot_main=try(gridExtra::grid.arrange(plot_io[["grid"]],plot_violin,plot_hist,plot_boxSideVar,plot_density,plot_capacity,
                           layout_matrix=rbind(c(1,2),
                                               c(3,3),
                                               c(4,5),
                                               c(6,6)
                           )),silent=FALSE)
    } else {
      plot_main=try(gridExtra::grid.arrange(plot_io[["grid"]],plot_violin,plot_hist,plot_boxSideVar,plot_density,plot_capacity,
                                        layout_matrix=rbind(c(1,2),
                                                            c(3,3),
                                                            c(4,5),
                                                            c(6,6),
                                                            c(6,6)
                                        )),silent=FALSE)
    }
  } else {
    if (is.null(cc_output$testing)){
      plot_main=try(gridExtra::grid.arrange(plot_io[["grid"]],plot_violin,plot_hist,plot_box,plot_density,plot_capacity,
                           layout_matrix=rbind(c(1,2),
                                               c(3,3),
                                               c(4,5),
                                               c(6,6)
                           )),silent=FALSE)
    } else {
      plot_main=try(gridExtra::grid.arrange(plot_io[["grid"]],plot_violin,plot_hist,plot_box,plot_density,plot_capacity,
                                        layout_matrix=rbind(c(1,2),
                                                            c(3,3),
                                                            c(4,5),
                                                            c(6,6),
                                                            c(6,6)
                                        )),silent=FALSE)
    }
  }
  
  if(!is.null(output_path)){
    try(ggplot2::ggsave(plot_main,file=paste(output_path,'MainPlot_full.pdf',sep=""),height=5*height,width=2*width,limitsize = FALSE),silent=FALSE)
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
      try(ggplot2::ggsave(plot_main_simp,file=paste(output_path,'MainPlot.pdf',sep=""),height=2*height,width=2*width,limitsize = FALSE),silent=FALSE)
    }
  } else {
    plot_main_simp=try(gridExtra::grid.arrange(plot_violin,plot_box,plot_capacity,
                                          layout_matrix=rbind(c(1,2),
                                                              c(3,3),
                                                              c(3,3),
                                                              c(3,3)
                                          )),silent=FALSE)
    if(!is.null(output_path)){
      try(ggplot2::ggsave(plot_main_simp,file=paste(output_path,'MainPlot.pdf',sep=""),height=3*height,width=2*width,limitsize = FALSE),silent=FALSE)
    }
  }
  
  
  
  
  graphOutput=list(plot_main,plot_io,plot_box,plot_violin,plot_hist,plot_boxSideVar,plot_capacity,plot_density,plot_main_simp)
}