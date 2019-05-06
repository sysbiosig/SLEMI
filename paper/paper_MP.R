### ### ### ### ### ### ### ### ### ### ### ### 
###                                         ###
### title: Script replicating results and   ### 
###        figures of MP of manuscript      ###
###                                         ###
### author: T. Jetka, T. Winarski,          ###
###         K. Nienaltowski, S.Blonski,     ###
###          M. Komorowski                  ###
### date: 06.05.2019                        ###
### Submission of manuscript:               ###
###                                         ###
### Information-theoretic analysis of       ###
### multivariate single-cell signaling      ###
### responses using SLEMI                   ###
###                                         ###
### ### ### ### ### ### ### ### ### ### ### ### 

#### Descritption ####
### ### ###
### The code below is divided into three sections (chunks; Use Rstudio's sections functionality for easy navigation):
### It should be run in the following order
### 1) Preliminary - setting up packages and working environment
### 2) Capacity - replicates Fig.1 A-B
### 3) Probabilities of discrimination - replicates Fig.1 C-D
### 
### Output will be saved in 'plotsMP/' under working directory (to be set in line 63)
##
## In default mode (5 repetition of bootstrap), running Capacity section takes approx. 2 hours with a single core.
## Set number of cores for parallel processing in line 83.
## Running Probability of discrimination section takes about 3 minutes.
##
## For a graphs like in MP (full diagnostic tests) set variable 'analysis_type' in line 84 to value 'full'
###
### Lines, where the activity of the user is needed, are emphasised with "FOR USER:"
### 
### Under the commands "print" we display the reference output. 
### Notice, that because of hardware differneces, sampling datasets 
### and probabilistic nature of optimisation procedures, results presented in script
### might be different. 
### However, it still enables to compare dimmensions and structure of the output.
###
### ### ###

#### Software configuration ####

# FOR USER: install packages that are required;
#           remember to load installed package
require("devtools")  # install.packages("devtools")
require("e1071")     # install.packages("e1071")
require("Hmisc")     # install.packages("Hmisc")
require("mvtnorm")   # install.packages("mvtnorm")
require("ggplot2")   # install.packages("ggplot2")
require("gridExtra") # install.packages("gridExtra")
require("corrplot") # install.packages("corrplot")


## 0. Install SLEMI package by: (if not installed already)
# install.packages("devtools") # if not installed
# library(devtools)
# install_github("sysbiosig/SLEMI")

#### Preliminaries ####
#1 Setting up working directory 
setwd("~/Downloads/") #FOR USER: Please change to a directory with downloaded scripts

#2. Loading neccessary libraries
source("aux_functions.R")
library(SLEMI)
aux_loadPackages(c("reshape2","ggplot2","stringr","doParallel","corrplot"))

#3 Inspect data saved in object data_nfkb, which stores our experimental measurements of NfkB signaling pathway
#Each row represents a single observation of a cell. Column 'signal' indicates the level of TNF stimulation for a given cell,
# while columns 'response_T', gives the normalised ratio of nuclear and cytoplasmic NfkB as described in SI. 
head(data_nfkb)

#4 Create output directory
data_type="Nfkb_5min_capacity"
full_path_out = paste('plotsMP/',data_type,'/',sep="")

#### Capacities ####
#5 Set parameters of algorithm - change for faster/more comprehensive analysis
random_seed=1234 # for reproducibility
set.seed(random_seed)
cores_num=1 # FOR USER: if possible, change to larger value to use parallel processing
analysis_type="short" 

bootstrap_prob=0.8
partition_trainfrac=0.6
if (analysis_type=="short"){
  bootstrap_num=10
  traintest_num=10
} else {
  bootstrap_num=50
  traintest_num=50
}


#6 Perform analysis of time-series data. It is a loop over all measured times: 0,3,6,...,120.
capacity_output<-list() #a list for storing capacity value
for (i_time in seq(from=0,to=120,by=3)) {
  capacity_chosen_times=seq(from=0,to=i_time,by=3)#
  capacity_ver=paste("capacity_TS_",i_time,sep="")
  path_capacity_output=paste(full_path_out,capacity_ver,"/",sep="")
  dir.create(path_capacity_output,recursive=TRUE)
  
  output_variables=paste("response",capacity_chosen_times,sep="_")
  Data_Analysis=data_nfkb[,c("signal",output_variables)]
  Data_Analysis$SignalNum=as.numeric(str_replace(str_match(as.character(Data_Analysis$signal),"[0-9]+[.]?[0-9]*"),",","."))
  Data_Analysis$SignalNum=as.factor(as.character(Data_Analysis[["SignalNum"]] ))
  
  Data_Analysis_final=Data_Analysis[ !apply(as.matrix(Data_Analysis[,c(output_variables) ]),1,function(x) any(is.na(x)|is.infinite(x) )), ]
  temp_capacity_output<-capacity_logreg_main(dataRaw=Data_Analysis_final,
                                             signal="SignalNum", response=output_variables, 
                                             output_path=path_capacity_output,data_out=FALSE,model_out=FALSE,
                                             testing=TRUE,scale=FALSE, 
                                             TestingSeed=random_seed,testing_cores=cores_num,
                                             boot_num=bootstrap_num,boot_prob=bootstrap_prob,
                                             traintest_num=traintest_num,partition_trainfrac=partition_trainfrac)
  
  capacity_output[[as.character(i_time)]]<- temp_capacity_output
}

saveRDS(capacity_output,file=paste(full_path_out,"capacity_TS.rds",sep="")) #saving capacity in a file

#7 Generate a graph for time series analysis
# It will be saved in "plotsMP/plot_1B.pdf" directory
dataTFT=capacity_output
dataPlot2=data.frame(capacity=sapply(dataTFT,function(x) mean(sapply(x$testing$bootstrap,function(y) y$cc)) ),
                     error1=sapply(dataTFT,function(x) quantile(sapply(x$testing$bootstrap,function(y) y$cc),probs=0.01)),
                     error2=sapply(dataTFT,function(x) quantile(sapply(x$testing$bootstrap,function(y) y$cc),probs=0.99)),
                     time=as.numeric(names(dataTFT)))
plot=ggplot(data=dataPlot2,aes(x=time,y=capacity))+
  geom_ribbon(aes(ymin=error1,ymax=error2),fill="grey75")+
  geom_line(colour="black",size=1.2)+
  theme_publ(version=2)+scale_y_continuous(limits=c(0,1.3),breaks = c(0,0.25,0.5,0.75,1,1.25))
filepath=paste("plotsMP/plot_1B.pdf",sep="")
ggsave(plot,file=filepath,height=6,width=8)



#8 Perform analysis of time-point data. It is a loop over all measured times: 0,3,6,...,120.
capacity_output=list() #list for storing capacity value
for (i_time in seq(from=0,to=120,by=3)) {
  capacity_chosen_times=i_time 
  capacity_ver=paste("capacity_TP_",i_time,sep="")
  path_capacity_output=paste(full_path_out,capacity_ver,"/",sep="")
  dir.create(path_capacity_output,recursive=TRUE)
  
  output_variables=paste("response",capacity_chosen_times,sep="_")
  Data_Analysis=data_nfkb[,c("signal",output_variables)]
  Data_Analysis$SignalNum=as.numeric(str_replace(str_match(as.character(Data_Analysis$signal),"[0-9]+[.]?[0-9]*"),",","."))
  Data_Analysis$SignalNum=as.factor(as.character(Data_Analysis[["SignalNum"]] ))
  
  Data_Analysis_final=Data_Analysis[ !apply(as.matrix(Data_Analysis[,c(output_variables) ]),1,function(x) any(is.na(x)|is.infinite(x) )), ]
  temp_capacity_output<-capacity_logreg_main(dataRaw=Data_Analysis_final,
                                             signal="SignalNum", response=output_variables, 
                                             output_path=path_capacity_output,data_out=FALSE,model_out=FALSE,
                                             testing=TRUE,scale=FALSE, 
                                             TestingSeed=random_seed,testing_cores=cores_num,
                                             boot_num=bootstrap_num,boot_prob=bootstrap_prob,
                                             traintest_num=traintest_num,partition_trainfrac=partition_trainfrac)
  
  capacity_output[[as.character(i_time)]]<- temp_capacity_output
}

saveRDS(capacity_output,file=paste(full_path_out,"capacityTP.rds",sep=""))  #saving capacity in a file


#9 Generate a graph for time-point analysis
# It will be saved in "plotsMP/plot_1A.pdf" directory
dataTSP=capacity_output
dataPlot2=data.frame(capacity=sapply(dataTSP,function(x) mean(sapply(x$testing$bootstrap,function(y) y$cc)) ),
                     error1=sapply(dataTSP,function(x) quantile(sapply(x$testing$bootstrap,function(y) y$cc),probs=0.01)),
                     error2=sapply(dataTSP,function(x) quantile(sapply(x$testing$bootstrap,function(y) y$cc),probs=0.99)),
                     time=as.numeric(names(dataTSP)))
plot=ggplot(data=dataPlot2,aes(x=time,y=capacity))+
  geom_ribbon(aes(ymin=error1,ymax=error2),fill = "grey75")+
  geom_line(colour="black",size=1.2)+
  theme_publ(version=2)+scale_y_continuous(limits=c(0,1.3),breaks = c(0,0.25,0.5,0.75,1,1.25))
filepath=paste("plotsMP/plot_1A.pdf",sep="")
ggsave(plot,file=filepath,height=6,width=8)



#### Probabilities of discrimination ####
# 10 Setting up enviornment
data_type="Nfkb_5min_probs_discrimination"
full_path_out = paste('plotsMP/',data_type,'/',sep="")

capacity_output<-list() #list for storing capacity value

chosen_stim=unique(data_nfkb$signal)
nstim=length(chosen_stim)

# 11 Estimate classificator
# loop for all input pairs
for (is in 1:(nstim-1) ){
  for (js in (is+1):nstim){
    capacity_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]]<-list()
    
    capacity_chosen_stim=chosen_stim[c(is,js)]
    signal_variable="signal"
    
    Data_Analysis=data_nfkb
    Data_Analysis=Data_Analysis[Data_Analysis[[signal_variable]] %in% capacity_chosen_stim,]
    Data_Analysis$SignalNum=as.numeric(str_replace(str_match(as.character(Data_Analysis[[signal_variable]]),"[0-9]+[.]?[0-9]*"),",","."))
    Data_Analysis$SignalNum=as.factor(as.character(Data_Analysis[["SignalNum"]] ))
    
    # time point
    capacity_chosen_times=21
    output_variable=paste("response",capacity_chosen_times,sep="_")

    capacity_ver=paste("Probs_",is,"_",js,"_TP",sep="")
    path_capacity_output=paste(full_path_out,capacity_ver,"/",sep="")
    dir.create(path_capacity_output,recursive=TRUE)
    
    Data_AnalysisFinal=Data_Analysis[ !apply(as.matrix(Data_Analysis[,c(output_variable) ]),1,function(x) any(is.na(x)|is.infinite(x) )), ]
    Data_Analysis_final=Data_AnalysisFinal[,c("SignalNum",output_variable)]
    
    capacity_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]][["tp"]]<-capacity_logreg_main(dataRaw=Data_Analysis_final,
                                                                                                     signal="SignalNum", response=output_variable,
                                                                                                     output_path=path_capacity_output,data_out=FALSE,
                                                                                                    testing=FALSE,scale=FALSE)
    
    #time series
    capacity_chosen_times=seq(from=0,to=120,by=3)
    output_variable=paste("response",capacity_chosen_times,sep="_")
    
    capacity_ver=paste("Probs_",is,"_",js,"_TS",sep="")
    path_capacity_output=paste(full_path_out,capacity_ver,"/",sep="")
    dir.create(path_capacity_output,recursive=TRUE)
    
    Data_AnalysisFinal=Data_Analysis[ !apply(as.matrix(Data_Analysis[,c(output_variable) ]),1,function(x) any(is.na(x)|is.infinite(x) )), ]
    Data_Analysis_final=Data_AnalysisFinal[,c("SignalNum",output_variable)]
    
    capacity_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]][["ts"]]<-capacity_logreg_main(dataRaw=Data_Analysis_final,
                                                                                                    signal="SignalNum", response=output_variable,
                                                                                                    output_path=path_capacity_output,data_out=FALSE,
                                                                                                    testing=FALSE,scale=FALSE)
  
    }
}
saveRDS(capacity_output,file = paste0(full_path_out,"capacity_probabilities.rds")) # save results in a file


#12 Extract results for time point data
ProbTP=matrix(0,nstim,nstim)

for (is in 1:(nstim-1) ){
  for (js in (is+1):nstim){
    w0_acc=c()
    temp_output=capacity_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]][["tp"]]
    temp_popt=temp_output$p_opt
    if (length(temp_popt)>2) {
      w0_acc=sum(temp_output$regression$overall[1])
    } else {
      w0_acc=sum(temp_output$regression$overall[1])
    }
    ProbTP[is,js]=w0_acc
   }
}


#13 Extract results for time series data
ProbTS=matrix(0,nstim,nstim)

for (is in 1:(nstim-1) ){
  for (js in (is+1):nstim){
    w0_acc=c()
    temp_output=capacity_output[[paste(chosen_stim[is],chosen_stim[js],sep="_")]][["ts"]]
    temp_popt=temp_output$p_opt
    if (length(temp_popt)>2) {
      w0_acc=sum(temp_output$regression$overall[1])
    } else {
      w0_acc=sum(temp_output$regression$overall[1])
    }
    ProbTS[is,js]=w0_acc
  }
}


#14 Visualise results as pie charts
#It will be saved in files "plotsMP/plot_1C.pdf" and "plotsMP/plot_1D.pdf"
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7","#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
ProbTP[ProbTP==0]=0.5
ProbTS[ProbTS==0]=0.5

pdf(paste0("plotsMP/plot_1C.pdf"),height=10,width=10)
corrplot(ProbTP,type = "upper", method = "pie",cl.lim = c(0.5, 1),col=col2(40))
dev.off()

pdf(paste0("plotsMP/plot_1D.pdf"),height=10,width=10)
corrplot(ProbTS,type = "upper", method = "pie",cl.lim = c(0.5, 1),col=col2(40))
dev.off()

#15 It can be also achieved without calculating capacities, by running function prob_discr_pairwise()
data_type="Nfkb_5min_probs_discrimination_2"
full_path_out = paste('plotsMP/',data_type,'/',sep="")

output_temp=prob_discr_pairwise(dataRaw=data_nfkb,signal="signal", response="response_21",output_path=paste0(full_path_out,"/time21/"))
output_temp=prob_discr_pairwise(dataRaw=data_nfkb,signal="signal", response=paste("response",seq(from=0,to=120,by=3),sep="_"),
                                output_path=paste0(full_path_out,"/time_TS/"))

