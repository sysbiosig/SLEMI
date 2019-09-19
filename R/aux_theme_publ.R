#' Tuning GGplot Themes
#'
#' Internal, auxillary functions
#'
#' This function changes the theme of plots created with the use of ggplot package
#' @param base_size integer that sets the default size of font used in the plot
#' @param base_family character that indicates what type of font should be used
#' @param version integer that changes the characteristic of the plot, values 1,2 and 3 accepted.
#' @keywords internal
#' @examples
#' library(ggplot2)
#' ggplot(data=data.frame(x=1:10,y=rnorm(10)),aes(x=x,y=y))+
#' geom_point()+SLEMI:::aux_theme_publ(version=2)
aux_theme_publ<-function (base_size = 12, base_family = "sans",version=1) { 
  if (version==2) {
    bgcolor = "default"
    bgcol <- ggthemes::ggthemes_data$hc$bg[bgcolor]
    ret <- ggplot2::theme(rect = ggplot2::element_rect(fill = bgcol, linetype = 0,colour = NA), 
                 text = ggplot2::element_text(size = base_size, family = base_family), 
                 title = ggplot2::element_text(size=20,hjust = 0.5,face="bold"), 
                 axis.title.x = ggplot2::element_text(size=18,hjust = 0.5,face="bold"), 
                 axis.title.y = ggplot2::element_text(size=18,hjust = 0.5,face="bold"),
                 axis.text = ggplot2::element_text(size=16,face="bold"),
                 axis.ticks.y=ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_line(),
                 axis.ticks.length=ggplot2::unit(0.2,"cm"),
                 axis.line = ggplot2::element_line(size = 1,colour="black",linetype="solid"),
                 panel.grid.major.y = ggplot2::element_line(color = "gray"), 
                 panel.grid.minor.y = ggplot2::element_blank(), 
                 panel.grid.major.x = ggplot2::element_blank(), 
                 panel.grid.minor.x = ggplot2::element_blank(), 
                 panel.border = ggplot2::element_blank(), 
                 panel.background = ggplot2::element_blank(), 
                 legend.position = "right")
  }
  ret
}