#' GGplot Themes
#' 
#' This function changes the theme of plots created with the use of ggplot package
#' @param base_size integer that sets the default size of font used in the plot
#' @param base_family character that indicates what type of font should be used
#' @param version integer that changes the characteristic of the plot, values 1,2 and 3 accepted.
#' @export
#' @examples
#' ggplot(data=data.frame(x=1:10,y=rnorm(10)),aes(x=x,y=y))+geom_point()+theme_publ(version=2)
theme_publ<-function (base_size = 12, base_family = "sans",version=1) { 
  if (version==1) {
  ret<-(ggthemes::theme_foundation(base_size = base_size, base_family = base_family) + 
          ggplot2::theme(line = ggplot2::element_line(), 
                rect = ggplot2::element_rect(fill = ggthemes::ggthemes_data$fivethirtyeight["ltgray"], linetype = 0, colour = NA), 
                text = ggplot2::element_text(colour = ggthemes::ggthemes_data$fivethirtyeight["dkgray"]), 
                axis.title = ggplot2::element_text(size=18,face="bold"),
                axis.title.y = ggplot2::element_text(angle=90),
                axis.text = ggplot2::element_text(size=14,face="bold"),
                axis.ticks = ggplot2::element_blank(),
                axis.line.x = ggplot2::element_line(),
                axis.line.y = ggplot2::element_blank(), 
                #legend.position="none",
                #panel.grid = element_line(colour = NULL), 
                panel.grid.major = ggplot2::element_line(colour = ggthemes::ggthemes_data$fivethirtyeight["medgray"]), 
                panel.grid.minor = ggplot2::element_blank(), 
                plot.title = ggplot2::element_text(hjust = 0, size = rel(1.75), face = "bold"), 
                plot.margin = ggplot2::unit(c(1,1, 1, 1), "lines"), 
                strip.background = ggplot2::element_rect(),
                strip.text = ggplot2::element_text(hjust = 0, size = rel(1), face = "bold") ))
} else if (version==2) {
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
} else if (version==3) {
  bgcolor = "default"
  bgcol <- ggthemes::ggthemes_data$hc$bg[bgcolor]
  ret <- ggplot2::theme(rect = ggplot2::element_rect(fill = bgcol, linetype = 0,colour = NA), 
               text = ggplot2::element_text(size = base_size, family = base_family), 
               title = ggplot2::element_text(size=20,hjust = 0.5,face="bold"), 
               axis.title.x = ggplot2::element_blank(),#element_text(size=18,hjust = 0.5,face="bold"), 
               axis.title.y = ggplot2::element_blank(),#element_text(size=18,hjust = 0.5,face="bold"),
               axis.text = ggplot2::element_text(size=16,face="bold"),
               axis.text.y = ggplot2::element_blank(),
               axis.ticks.y=ggplot2::element_blank(),
               axis.ticks.x = ggplot2::element_line(),
               axis.ticks.length=ggplot2::unit(0.2,"cm"),
               axis.line = ggplot2::element_line(size = 1,colour="black",linetype="solid"),
               axis.line.y = ggplot2::element_blank(),
               panel.grid.major.y = ggplot2::element_blank(), 
               panel.grid.minor.y = ggplot2::element_blank(), 
               panel.grid.major.x = ggplot2::element_blank(), 
               panel.grid.minor.x = ggplot2::element_blank(), 
               panel.border = ggplot2::element_blank(), 
               panel.background = ggplot2::element_blank(), 
               legend.position = "right")
}  else {
  bgcolor = "darkunica"
  bgcol <- ggthemes::ggthemes_data$hc$bg[bgcolor]
  ret <- ggplot2::theme(rect = ggplot2::element_rect(fill = bgcol, linetype = 0,colour = NA), 
               text = ggplot2::element_text(size = base_size, family = base_family), 
               title = ggplot2::element_text(size=20,hjust = 0.5,face="bold"), 
               axis.title.x = ggplot2::element_text(size=18,hjust = 0.5,face="bold"), 
               axis.title.y = ggplot2::element_text(size=18,hjust = 0.5,face="bold"),
               axis.text = ggplot2::element_text(size=16,face="bold"),
               panel.grid.major.y = ggplot2::element_line(color = "gray"), 
               panel.grid.minor.y = ggplot2::element_blank(), 
               panel.grid.major.x = ggplot2::element_blank(), 
               panel.grid.minor.x = ggplot2::element_blank(), 
               panel.border = ggplot2::element_blank(), 
               panel.background = ggplot2::element_blank(), 
               legend.position = "none")
  ret <- (ret + ggplot2::theme(rect = ggplot2::element_rect(fill = bgcol), 
                      text = ggplot2::element_text(colour = "#A0A0A3"), 
                      title = ggplot2::element_text(colour = "#FFFFFF"), 
                      axis.title.x = ggplot2::element_text(colour = "#A0A0A3"), 
                      axis.title.y = ggplot2::element_text(colour = "#A0A0A3"), 
                      panel.grid.major.y = ggplot2::element_line(color = "gray"), 
                      legend.title = ggplot2::element_text(colour = "#A0A0A3")))
}
  ret
}