#' Plotting output of capacity estimation. Auxillary functions.
#' 
#' INPUT:
#' @param data is a data.frame object
#' @param signal is a character object that indicates columns of data that should be treated as channel's input
#' @param response is a character vector that indicates columns of data that should be treated as channel's output
#' @param side_variables is a character vector that indicates side variables' columns of data
#' @param cc_output a list that is the output of capacity_logreg_algorithm function
#' @param path character giving the directory, where graphs should be saved
#' @param height integer indicating the height of a single plot
#' @param width integer indicating the width of a single plot
#' @keywords internal
#' @examples 
#' 
capacity_output_graph_io<-function(data,signal,response,path,height,width){
  
  data_colnames=colnames(data)
  response_length=length(response)
  signalNUM=paste(signal,"_RAW",sep="")
  
  if (any(data_colnames %in% (signalNUM))){

    maxSignal=max(data[[signalNUM]])
    
    dataPlot=reshape2::melt(data[,c(signalNUM,response)],id.vars=c(signalNUM))
    plot<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x=signalNUM,y="value"))+ggplot2::stat_summary(fun.y=mean,geom="line",size=1.35)+
      ggplot2::facet_grid(.~variable)+
      ggplot2::scale_y_continuous(paste("Response ",sep="") )+ggplot2::scale_x_continuous(paste("Signal - ",signal,sep=""),limits=c(0,maxSignal*1.05) )+
      ggplot2::ggtitle("Input-Output Relation")+
      theme_publ(version=2)
    
    plot_col<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x=signalNUM,y="value"),ggplot2::aes(colour=factor(variable)))+ggplot2::stat_summary(fun.y=mean,geom="line",size=1.35)+
      ggplot2::scale_y_continuous(paste("Response ",sep="") )+ggplot2::scale_x_continuous(paste("Signal - ",signal,sep=""),limits=c(0,maxSignal*1.05) )+
      ggplot2::ggtitle("Input-Output Relation")+
      theme_publ(version=2)
  } else {
    dataPlot=reshape2::melt(data[,c(signal,response)],id.vars=c(signal))
    plot<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x=signal,y="value"))+
      ggplot2::stat_summary(fun.y=mean,geom="point",size=2.35)+
      ggplot2::facet_grid(.~variable)+
      ggplot2::scale_y_continuous(paste("Response ",sep="") )+
      ggplot2::scale_x_discrete(paste("Signal - ",signal,sep="") )+
      ggplot2::ggtitle("Input-Output Relation")+
      theme_publ(version=2)
    
    plot_col<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x=signal,y="value"),ggplot2::aes(colour=factor(variable)))+
      ggplot2::stat_summary(fun.y=mean,geom="point",size=2.35)+
      ggplot2::scale_y_continuous(paste("Response ",sep="") )+
      ggplot2::scale_x_discrete(paste("Signal - ",signal,sep="") )+
      ggplot2::ggtitle("Input-Output Relation")+
      theme_publ(version=2)
  }
  
  ggplot2::ggsave(plot,file=paste(path,'io_relation.pdf',sep=""),height=height,width=min(49.5,response_length*width))
  #ggplot2::ggsave(plot_col,file=paste(path,'io2_relation.pdf',sep=""),height=height,width=width)
  
  list(grid=plot,colour=plot_col)
}


#' @rdname capacity_output_graph_io
capacity_output_graph_boxplots<-function(data,signal,response,path,height=4,width=6){
  
  data_colnames=colnames(data)
  response_length=length(response)
  signalNUM=paste(signal,"_RAW",sep="")
  
  if (any(data_colnames %in% (signalNUM))){
  signal=signalNUM
    }
  
  if (is.numeric(data[,signal])){
    data[,signal]<-factor(data[,signal])
  }
  
  dataPlot=reshape2::melt(data[,c(signal,response)],id.vars=c(signal))
  plot<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x=signal,y="value") ) + 
    ggplot2::geom_boxplot()+# geom_point() +
    ggplot2::facet_grid(.~variable)+
    ggplot2::scale_y_continuous(name=paste("Output",sep="") )+
    #ggplot2::scale_x_discrete(name=paste("Input",signal,sep="") )+
    ggplot2::scale_x_discrete(name=paste("Input",sep="") )+
    ggplot2::ggtitle("Boxplots")+
    theme_publ(version=2)
  
  ggplot2::ggsave(plot,file=paste(path,'boxplots.pdf',sep=""),height=height,width=width)
  
  plot
}


#' @rdname capacity_output_graph_io
capacity_output_graph_violinMean<-function(data,signal,response,path,height=4,width=6){
  
  data_colnames=colnames(data)
  response_length=length(response)
  signalNUM=paste(signal,"_RAW",sep="")
  
  if (any(data_colnames %in% (signalNUM))){
    
    maxSignal=max(data[[signalNUM]])
    
    dataPlot=reshape2::melt(data[,c(signalNUM,response)],id.vars=c(signalNUM))
    plot<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x=signalNUM,y="value"))+ggplot2::geom_violin(ggplot2::aes_string(group=signalNUM),scale="width")+
      ggplot2::stat_summary(fun.y=mean,geom="line",size=1.15)+ggplot2::stat_summary(fun.y=mean,geom="point",size=1.5)+
      ggplot2::facet_grid(variable~.)+
      ggplot2::scale_y_continuous(paste("Output",sep="") )+
      #ggplot2::scale_x_continuous(paste("Input",signal,sep=""),limits=c(0,maxSignal*1.05) )+
      ggplot2::scale_x_continuous(paste("Input",sep=""),limits=c(0,maxSignal*1.05) )+
      ggplot2::ggtitle("ViolinPlot")+
      theme_publ(version=2)
  } else {
    dataPlot=reshape2::melt(data[,c(signal,response)],id.vars=c(signal))
    plot<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x=signal,y="value"))+ggplot2::geom_violin(ggplot2::aes_string(group=signal),scale="width")+
      ggplot2::stat_summary(fun.y=mean,geom="point",size=1.5)+
      ggplot2::facet_grid(variable~.)+
      ggplot2::scale_y_continuous(paste("Output",sep="") )+
      #ggplot2::scale_x_discrete(paste("Input",signal,sep="") )+
      ggplot2::scale_x_discrete(paste("Input",sep="") )+
      ggplot2::ggtitle("ViolinPlot")+
      theme_publ(version=2)
  }
  
  ggplot2:: ggsave(plot,file=paste(path,'violin.pdf',sep=""),height=height,width=min(49.5,response_length*width))
  
  plot
}


#' @rdname capacity_output_graph_io
capacity_output_graph_histograms<-function(data,signal,response,path,height=4,width=6){
  
  response_length=length(response)
  
  obsv_num=unlist(by(data.frame(data[,response]),data[[signal]],function(x) nrow(x)))
  minrange = min(unlist(sapply(by(data.frame(data[,response]),data[[signal]],function(x) sapply(x,range) ),function(x) x[2,]-x[1,])))
  
  tempformula=as.formula(paste("variable ~ ",signal,sep=""))
  
  dataPlot=reshape2::melt(data[,c(signal,response)],id.vars=c(signal))
  plotTemp<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x="value"),ggplot2::aes(y=..density..) ) + 
    ggplot2::geom_histogram(binwidth=minrange/10,ggplot2::aes(y=..density..)) + 
    ggplot2::facet_grid(tempformula)+
    ggplot2::coord_flip()
  
  plotText0=grid::textGrob("Number of observations:")
  plotText=grid::textGrob(paste(c(paste(names(obsv_num),":",obsv_num,sep="")),collapse="   "))
  plot=gridExtra::grid.arrange(plotText0,plotText,plotTemp,layout_matrix=matrix(c(1,2,3,3,3,3),nrow=6,ncol=1))
  ggplot2::ggsave(plot,file=paste(path,'histograms.pdf',sep=""),height=min(49.5,response_length*height),width=min(49.5,3*width) )
  
  plot  
}
# capacity_output_graph_histograms<-function(data,signal,response,path,height=4,width=6){
#   response_length=length(response)
#   obsv_num=unlist(by(data.frame(data[,response]),data[[signal]],function(x) nrow(x)))
#   minrange = min(unlist(sapply(by(data.frame(data[,response]),data[[signal]],function(x) sapply(x,range) ),function(x) x[2,]-x[1,])))
#   tempformula=as.formula(paste("variable ~ ",signal,sep=""))
#   dataPlot=reshape2::melt(data[,c(signal,response)],id.vars=c(signal))
#   plotTemp<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x="value"),ggplot2::aes(y=..density..) ) + 
#     ggplot2::geom_histogram(binwidth=minrange/10,ggplot2::aes(y=..density..)) + 
#     ggplot2::facet_grid(tempformula)+
#     ggplot2::coord_flip()
#   r<-print(plotTemp)
#   dc <- plyr::dlply(r$data[[1]], plyr::.(PANEL), function(x) max(x$density))
#   tmp_data=data.frame(y=max(unlist(dc))/2,x=Inf,lab= paste("No. obs:", obsv_num,sep=" "), signal=names(obsv_num) )
#   colnames(tmp_data)[4]<-signal
#   plot<-plotTemp + ggplot2::geom_text(data=tmp_data,ggplot2::aes(x=x,y=y,label=lab),vjust=1)+    
#     ggplot2::scale_y_continuous(paste("Output",sep="") )+
#     ggplot2::ggtitle("Histograms")+
#     theme_publ(version=2)
#   ggplot2::ggsave(plot,file=paste(path,'histograms.pdf',sep=""),height=min(49.5,response_length*height),width=min(49.5,3*width) )
#   plot  
# }


#' @rdname capacity_output_graph_io
capacity_output_graph_boxplotsSideVar<-function(data,signal,side_variables,path,height=4,width=6){
  
  if (is.null(side_variables)) {
    plot=grid::textGrob(" ")  
  } else if (!all(sapply(data[,side_variables],function(x) is.numeric(x) ))) {
    plot=grid::textGrob(" ")  
  } else  {
    dataTemp=data[,colnames(data)%in% c(signal,side_variables)]
    dataMelt=reshape2::melt(dataTemp, id.vars=c(signal))
    
    plot<-ggplot2::ggplot(data=dataMelt,ggplot2::aes_string(x=signal,y="value") ) + ggplot2::geom_boxplot()+ ggplot2::facet_grid(variable~.) +
      ggplot2::scale_y_continuous("side variable")+
      ggplot2::scale_x_discrete(paste("Input",signal,sep="") )+
      ggplot2::ggtitle("Boxplots")+
      theme_publ(version=2)
    
    ggplot2::ggsave(plot,file=paste(path,'boxplots_sideVar.pdf',sep=""),height=height,width=width)
  }
  
  plot
}


#' @rdname capacity_output_graph_io
capacity_output_graph_capacity<-function(cc_output,path,height=4,width=6){
  
  if (is.null(cc_output$testing)){
    
    plot<-ggplot2::ggplot(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y) ) + 
      ggplot2::geom_point(size=4,shape=15) + 
      ggplot2::geom_text(ggplot2::aes(label=round(x,digits=3),y=0.075),size=8)+
      ggplot2::geom_line(data=data.frame(x=seq(from=0.9*cc_output$cc,to=1.1*cc_output$cc,length=10),y=0),arrow = ggplot2::arrow())+
      ggplot2::scale_y_continuous("",limits=c(-0.05,0.1))+ggplot2::scale_x_continuous( "Capacity (bits)" )+
      ggplot2::ggtitle("Capacity")+
      theme_publ(version=2)+
      ggplot2::theme(axis.ticks.y=ggplot2::element_blank(),
            axis.text.y=ggplot2::element_blank(),
            axis.title.y=ggplot2::element_blank(),
            panel.grid.major.y=ggplot2::element_blank(),
            axis.line.y=ggplot2::element_blank())
    
    ggplot2::ggsave(plot,file=paste(path,'capacity.pdf',sep=""),height=2,width=6)
    
  } else if (length(cc_output$testing)==2) {
    
    boot_results=sapply(cc_output$testing$bootstrap,function(x) x$cc)
    boot_sd=sd(boot_results)
    boot_mean=mean(boot_results)
    
    plot1<-ggplot2::ggplot(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y) ) + 
      ggplot2::geom_point(size=6,shape=15) + 
      ggplot2::geom_text(ggplot2::aes(label=round(x,digits=3),y=0.075),size=8)+
      ggplot2::geom_line(data=data.frame(x=seq(from=0.95*(cc_output$cc-boot_sd),to=1.05*(cc_output$cc+boot_sd),length=10),y=0),arrow = ggplot2::arrow())+
      ggplot2::geom_point(data=data.frame(x=c(boot_mean-boot_sd,boot_mean+boot_sd),y=c(0,0)),shape=124,colour="red",size=10)+ 
      ggplot2::geom_text(data=data.frame(x=c(0.985*(boot_mean-boot_sd),1.015*(boot_mean+boot_sd)),y=c(0.05,0.05)),ggplot2::aes(label=round(x,digits=3),y=0.025),size=6, colour="red")+
      ggplot2::geom_point(data=data.frame(x=c(boot_mean),y=c(0)),shape=15,colour="red",size=4)+
      ggplot2::scale_y_continuous("",limits=c(-0.05,0.1))+
      ggplot2::scale_x_continuous( "Capacity (bits)",breaks=c(0.95*(cc_output$cc-boot_sd),1.05*(cc_output$cc+boot_sd)) )+
      ggplot2::ggtitle("Capacity")+
      theme_publ(version=2)+
      ggplot2::theme(axis.ticks.y=ggplot2::element_blank(),
            axis.text.y=ggplot2::element_blank(),
            axis.title.y=ggplot2::element_blank(),
            panel.grid.major.y=ggplot2::element_blank(),
            axis.line.y=ggplot2::element_blank())
    
    plot2<-ggplot2::ggplot( data=data.frame(x=sapply(cc_output$testing$bootstrap,function(x) x$cc)),ggplot2::aes(x=x,y=..density..) ) + 
      ggplot2::geom_histogram( )+
      ggplot2::geom_point(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y),size=4,shape=15) +
      ggplot2::ggtitle("Bootstrap")+ggplot2::scale_x_continuous( "Capacity (bits)" )+
      theme_publ(version=2)
    
    plot4<-ggplot2::ggplot( data=data.frame(x=sapply(cc_output$testing$traintest,function(x) x$cc)),ggplot2::aes(x=x,y=..density..) ) + 
      ggplot2::geom_histogram( )+
      ggplot2::geom_point(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y),size=4,shape=15) +
      ggplot2::ggtitle("TrainTest")+ggplot2::scale_x_continuous( "Capacity (bits)" )+
      theme_publ(version=2)
    
    plot5=grid::textGrob(paste("PV-left:",round(cc_output$testing_pv$bootstrap[1],digits=3),sep=" " ))
     plot6=grid::textGrob(paste("PV-right:",round(cc_output$testing_pv$bootstrap[2],digits=3),sep=" " )) 

    plot9=grid::textGrob(paste("PV-left:",round(cc_output$testing_pv$traintest[1],digits=3),sep=" " ))
    plot10=grid::textGrob(paste("PV-right:",round(cc_output$testing_pv$traintest[2],digits=3),sep=" " )) 

    
    plot=gridExtra::grid.arrange(plot1,plot2,plot4,plot5,plot6,plot9,plot10,
                      layout_matrix=rbind(c(1,1),c(1,1),c(2,4),c(2,5),c(3,6),c(3,7)))
    
    
    ggplot2::ggsave(plot,file=paste(path,'capacity.pdf',sep=""),height=1.5*height,width=width)
    
  } else {
    boot_results=sapply(cc_output$testing$bootstrap,function(x) x$cc)
    boot_sd=sd(boot_results)
    boot_mean=mean(boot_results)
    
    plot1<-ggplot2::ggplot(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y) ) + 
      ggplot2::geom_point(size=6,shape=15) + 
      ggplot2::geom_text(ggplot2::aes(label=round(x,digits=3),y=0.075),size=8)+
      ggplot2::geom_line(data=data.frame(x=seq(from=0.95*(cc_output$cc-boot_sd),to=1.05*(cc_output$cc+boot_sd),length=10),y=0),arrow = ggplot2::arrow())+
      ggplot2::geom_point(data=data.frame(x=c(boot_mean-boot_sd,boot_mean+boot_sd),y=c(0,0)),shape=124,colour="red",size=10)+ 
      ggplot2::geom_text(data=data.frame(x=c(0.985*(boot_mean-boot_sd),1.015*(boot_mean+boot_sd)),y=c(0.05,0.05)),ggplot2::aes(label=round(x,digits=3),y=0.025),size=6, colour="red")+
      ggplot2::geom_point(data=data.frame(x=c(boot_mean),y=c(0)),shape=15,colour="red",size=4)+
      ggplot2::scale_y_continuous("",limits=c(-0.05,0.1))+
      ggplot2::scale_x_continuous( "Capacity (bits)",breaks=c(0.95*(cc_output$cc-boot_sd),1.05*(cc_output$cc+boot_sd)) )+
      ggplot2::ggtitle("Capacity")+
      theme_publ(version=2)+
      ggplot2::theme(axis.ticks.y=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     panel.grid.major.y=ggplot2::element_blank(),
                     axis.line.y=ggplot2::element_blank())
    
    plot2<-ggplot2::ggplot( data=data.frame(x=sapply(cc_output$testing$bootstrap,function(x) x$cc)),ggplot2::aes(x=x,y=..density..) ) + 
      ggplot2::geom_histogram( )+
      ggplot2::geom_point(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y),size=4,shape=15) +
      ggplot2::ggtitle("Bootstrap")+ggplot2::scale_x_continuous( "Capacity (bits)" )+
      theme_publ(version=2)
    
    plot3<-ggplot2::ggplot( data=data.frame(x=sapply(cc_output$testing$resamplingMorph,function(x) x$cc)),ggplot2::aes(x=x,y=..density..) ) + 
      ggplot2::geom_histogram( )+
      ggplot2::geom_point(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y),size=4,shape=15) +
      ggplot2::ggtitle("SideVar Resampling")+ggplot2::scale_x_continuous( "Capacity (bits)" )+
      theme_publ(version=2)
    
    plot4<-ggplot2::ggplot( data=data.frame(x=sapply(cc_output$testing$traintest,function(x) x$cc)),ggplot2::aes(x=x,y=..density..) ) + 
      ggplot2::geom_histogram( )+
      ggplot2::geom_point(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y),size=4,shape=15) +
      ggplot2::ggtitle("TrainTest")+ggplot2::scale_x_continuous( "Capacity (bits)" )+
      theme_publ(version=2)
    
    plot3b<-ggplot2::ggplot( data=data.frame(x=sapply(cc_output$testing$bootResampMorph,function(x) x$cc)),ggplot2::aes(x=x,y=..density..) ) + 
      ggplot2::geom_histogram( )+
      ggplot2::geom_point(data=data.frame(x=cc_output$cc,y=0),ggplot2::aes(x=x,y=y),size=4,shape=15) +
      ggplot2::ggtitle("Boot SideVar Resampling")+ggplot2::scale_x_continuous( "Capacity (bits)" )+
      theme_publ(version=2)
    
    plot5=grid::textGrob(paste("PV-left:",round(cc_output$testing_pv$bootstrap[1],digits=3),sep=" " ))
    plot6=grid::textGrob(paste("PV-right:",round(cc_output$testing_pv$bootstrap[2],digits=3),sep=" " )) 
    plot7=grid::textGrob(paste("PV-left:",round(cc_output$testing_pv$resamplingMorph[1],digits=3),sep=" " ))
    plot8=grid::textGrob(paste("PV-right:",round(cc_output$testing_pv$resamplingMorph[2],digits=3),sep=" " )) 
    plot9=grid::textGrob(paste("PV-left:",round(cc_output$testing_pv$traintest[1],digits=3),sep=" " ))
    plot10=grid::textGrob(paste("PV-right:",round(cc_output$testing_pv$traintest[2],digits=3),sep=" " )) 
    plot11=grid::textGrob(paste("PV-left:",round(cc_output$testing_pv$bootResampMorph[1],digits=3),sep=" " ))
    plot12=grid::textGrob(paste("PV-right:",round(cc_output$testing_pv$bootResampMorph[2],digits=3),sep=" " )) 
    
    plot=gridExtra::grid.arrange(plot1,plot2,plot3,plot4,plot3b,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12,
                                 layout_matrix=rbind(c(1,1),c(1,1),c(2,6),c(2,7),c(3,8),c(3,9),c(4,10),c(4,11),c(5,12),c(5,13)))
    
    
    ggplot2::ggsave(plot,file=paste(path,'capacity.pdf',sep=""),height=2*height,width=width)
  }
  
  plot
}



#' @rdname capacity_output_graph_io
capacity_output_graph_densities<-function(data,signal,response,path,height=4,width=6){
  
  dataPlot=reshape2::melt(data[,c(signal,response)],id.vars=c(signal))
  plot<-ggplot2::ggplot(data=dataPlot,ggplot2::aes_string(x="value",fill=signal,colour=signal)) + ggplot2::geom_density(alpha = 0.3)+
    ggplot2::facet_grid(variable~.)+
    ggplot2::scale_x_continuous(paste("Output",sep="") )+ ggplot2::scale_y_continuous("Density")+
    ggplot2::ggtitle("Density")+
    theme_publ(version=2)+ggplot2::theme(legend.position = "right")
  
  ggplot2::ggsave(plot,file=paste(path,'kdensities.pdf',sep=""),height=height,width=width)
  
  plot
}
