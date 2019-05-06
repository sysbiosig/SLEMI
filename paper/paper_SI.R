### ### ### ### ### ### ### ### ### ### ### ### 
###                                         ###
### title: Script replicating results and   ### 
###        figures of SI of manuscript      ###
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
### 2) Comparison - replicates Fig. S1 - shows the comaprison of our method to the KNN approach.
### 3) Validation - replicates Fig. S2 - shows the performance of our method in four examples of simple channels
## In default mode (10 repetition of data sampling), running Validation section takes 1 hour, 
## similarly computations in Comparison section also take approximately 1 hour.
## with a single core. Set number of cores for parallel processing in line 76.
##
## For a graphs like in SI (full diagnostic tests) set variable 'analysis_type' in line 77 to value 'full'
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
# remember to load installed package
require("devtools")  # install.packages("devtools")
require("e1071")     # install.packages("e1071")
require("mvtnorm")   # install.packages("mvtnorm")
require("ggplot2")   # install.packages("ggplot2")
require("gridExtra") # install.packages("gridExtra")

# for knn method
require("nloptr")   # install.packages("nloptr")
require("FNN")      # install.packages("FNN")
require("DEoptimR") # install.packages("DEoptimR")
require("TDA")      # install.packages("TDA")


## 0. Install SLEMI package by:
# install.packages("devtools") # if not installed
# library(devtools)
# install_github("sysbiosig/SLEMI")

#### 1. Preliminaries ####
#1.1 Setting up working directory 
setwd("~/Downloads/") #FOR USER: Please change to a directory with downloaded scipts

#1.2 Loading neccessary libraries
source("aux_functions.R")
library(SLEMI)
aux_loadPackages(c("reshape2","ggplot2","stringr","doParallel","gridExtra"))

#1.3 Create output directory
dir.create("plotsSI/synthetic_validation/",recursive = TRUE)
dir.create("plotsSI/synthetic_comparison/",recursive = TRUE)

#1.4 Set general parameters
boot_cores=3 # FOR USER: if possible, change to larger value to use parallel processing
analysis_type="short" 
#analysis_type="long" #FOR USER: uncomment this line for an exact replication of the results
random_seed=12345  # for reproducibility
set.seed(random_seed)

if (analysis_type=="short"){
  boot_num=3
} else {
  boot_num=40
}


#### 2. Comparison ####
## After running this chunk, plots similar to those in Fig. S2 will be created in directory plotsSI/synthetic_comparison/
## Without any changes, it will take ~1hour to re-calculate all examples
## 2.1 Loading additional packages (for using KNN method)
aux_loadPackages(c("nloptr","FNN","DEoptimR","TDA"))

##Influence of dimension of the output ####

## 2.2 Seting parameters of estimation
ns=c(500,2000,4000)
sds=c(1)
corr_init=0
if (analysis_type=="short"){
  i_outputs=paste(seq(from=2,to=30,by=7),"d",sep="")
} else {
  i_outputs=paste(2:30,"d",sep="")
} 
i_class="2"
i_type="multigauss_dimensions"

path_output_main=paste('plotsSI/synthetic_comparison/',i_type,"/",sep="")
dir.create(path_output_main,recursive = TRUE)

OutputList=list()
OutputListT=list()

# 2.3 Setting specifics of the channel
classes_num = as.numeric(i_class) 

i_output="1d"
dimension_num=as.numeric(str_match(i_output,"[0-9]+"))

mean_diff_init=c(2,rep(0,(dimension_num-1) ))
mean_init=rep(0,dimension_num)
sd_init=rep(1,dimension_num)

example_means1 = mean_init
example_means2 = mean_init+mean_diff_init

example_sigma1 = diag(sd_init)
example_sigma2 = diag(sd_init)

# 2.4 Finding exact solution
tempoutput_exact <- cc_exact(type="gauss",classes=i_class,
                             dimension= dimension_num,
                             par=list(mu1=example_means1,
                                      mu2=example_means2,
                                      sd1=example_sigma1,
                                      sd2=example_sigma2),
                             x_length=10000,s_length=200)


# 2.5 Estimating capacity by LogReg and KNN methods for different dimensions of the output
for (i_output in i_outputs){
  
  classes_num = as.numeric(i_class) 
  dimension_num=as.numeric(str_match(i_output,"[0-9]+"))
  
  mean_diff_init=c(2,rep(0,(dimension_num-1) ))
  mean_init=rep(0,dimension_num)
  sd_init=rep(1,dimension_num)
  
  example_means1 = mean_init
  example_means2 = mean_init+mean_diff_init
  
  example_sigma1 = diag(sd_init)
  example_sigma2 = diag(sd_init)
  
  OutputList[[as.character(i_output)]]=list()
  OutputList[[as.character(i_output)]]$exact=tempoutput_exact
  OutputListT[[as.character(i_output)]]=list()
  
  for (i_n in ns){
    print(paste(i_n,i_output,sep="_"))
    path_output_single=paste(path_output_main,
                             'dim_',dimension_num,'_N_',i_n,'/',
                             sep="")
    
    OutputList[[as.character(i_output)]][[as.character(i_n)]]=list()  
    OutputListT[[as.character(i_output)]][[as.character(i_n)]]=list()
    
    
    ##SLEMI algorithm
    t1=proc.time()
    cl=parallel::makeCluster(boot_cores)
    doParallel::registerDoParallel(cl)
    OutList=foreach::foreach(i=1:boot_num,.packages = c("SLEMI","mvtnorm","FNN","TDA","DEoptimR")) %do% {
      
      tempdata = data.frame(signal = c(t(replicate(i_n,LETTERS[1:classes_num]))) ,
                            rbind(
                              rmvnorm(i_n,mean=example_means1,sigma=example_sigma1),
                              rmvnorm(i_n,mean=example_means2,sigma=example_sigma2) 
                            ))
      
      path_output_single_boot=paste(path_output_single,i,"/",sep="")
      dir.create(path_output_single_boot,recursive = TRUE)
      tempoutput  <- capacity_logreg_main(dataRaw=tempdata,
                                          signal="signal", response=paste("X",1:dimension_num,sep=""),
                                          formula_string = paste("signal~", paste(paste("X",1:dimension_num,sep=""),collapse="+"),sep=""),
                                          output_path=path_output_single_boot,plot_height=8,plot_width=12,
                                          testing=FALSE,model_out = FALSE,data_out = FALSE)
      tempout=tempoutput$cc
      rm(tempdata,tempoutput)
      tempout
    }
    stopCluster(cl)
    
    OutputListT[[as.character(i_output)]][[as.character(i_n)]]$lr= (proc.time()-t1)[3]
    OutputList[[as.character(i_output)]][[as.character(i_n)]]$lr=OutList 
    
    ##KNN algorithm
    t1=proc.time()
    cl=parallel::makeCluster(boot_cores)
    doParallel::registerDoParallel(cl)
    OutList=foreach::foreach(i=1:boot_num,.packages = c("SLEMI","mvtnorm","FNN","TDA","DEoptimR")) %do% {
      
      tempdata = data.frame(signal = c(t(replicate(i_n,LETTERS[1:classes_num]))) ,
                            rbind(
                              rmvnorm(i_n,mean=example_means1,sigma=example_sigma1),
                              rmvnorm(i_n,mean=example_means2,sigma=example_sigma2) 
                            ))
      
      ni_class=as.numeric(i_class)
      
      knndens=list()
      for (is in LETTERS[1:classes_num]){
        knndens[[is]]=knnDE(tempdata[tempdata$signal==is,-1],tempdata[,-1],10)
      }
      
      eval_f_args<-function(q,data,knndens){
        
        signals=as.character(unique(data[,1]))
        fri=lapply(signals,function(x) { knndens[[x]][data[,1]==x]  } )
        fr=Reduce("+",mapply(function(x,y){ y*knndens[[x]] },signals,q,SIMPLIFY=FALSE))
        
        hrs=-sum( mapply(function(x,y){
          sum(log2(x))*y*(1/length(x))
        },fri,q,SIMPLIFY=TRUE)  )
        
        hr=-sum( mapply(function(x,y){
          sum(log2(fr[data[1]==x]))*y*(1/length(x[data[1]==x]))
        },signals,q,SIMPLIFY=TRUE)  )
        
        return(hrs-hr)
      }
      
      
      eval_f<-function(q){return(eval_f_args(q,data=tempdata,knndens= knndens))}
      eval_g<-function(q){return(1-sum(q))}
      
      resDE<-JDEoptim(rep(0,ni_class),rep(1,ni_class),
                      fn = eval_f, constr = eval_g, meq = 1,
                      tol = 1e-6, trace = FALSE, triter = 50,maxiter = 6000)
      
      tempout=-resDE$value
      
      rm(tempdata,knndens)
      tempout
    }
    stopCluster(cl)
    
    OutputListT[[as.character(i_output)]][[as.character(i_n)]]$knn= (proc.time()-t1)[3]
    OutputList[[as.character(i_output)]][[as.character(i_n)]]$knn=OutList 
    
    
  }
}

saveRDS(OutputListT,file=paste0(path_output_main,"/times.rds")) #saving objects in files
saveRDS(OutputList,file=paste0(path_output_main,"/knn.rds"))


# 2.6 Generating graphs
dir.create("plotsSI/synthetic_comparison/",recursive = TRUE)
df=melt(lapply(OutputList,function(x) lapply(x[c("500","2000","4000")],function(y){ lapply(y,function(z) z) })))
colnames(df)<-c("capacity","boot","method","sampleSize","dimension")
df$dimension=factor(df$dimension,levels=i_outputs)
df$sampleSize=factor(df$sampleSize,levels=c("500","2000","4000"))

plot=ggplot(data=df,aes(x=dimension,y=capacity,colour=method))+facet_grid(sampleSize~.)+
  stat_summary(fun.data="mean_sdl")+
  theme_publ(version=2)+geom_hline(yintercept = OutputList$`2d`$exact$cc)
ggsave(plot,file="plotsSI/plot_S1A.pdf",height=6,width=8)


dft=melt(lapply(OutputListT,function(x) lapply(x[c("500","2000","4000")],function(y){ lapply(y,function(z) z) })))
colnames(dft)<-c("computationTime","method","sampleSize","dimension")
dft$dimension=factor(dft$dimension,levels=i_outputs)
dft$dimensionNum=as.numeric(str_match(dft$dimension,"[0-9]+"))
dft$sampleSize=factor(dft$sampleSize,levels=c("500","2000","4000"))
dft$computationTime=dft$computationTime

plot=ggplot(data=dft,aes(x=dimensionNum,y=computationTime,colour=method))+
  facet_grid(sampleSize~.)+geom_line(size=2)+
  scale_x_continuous("dimension")+scale_y_continuous("capacity")+
  theme_publ(version=2)
ggsave(plot,file="plotsSI/plot_S1C.pdf",height=6,width=8)




## 2.7 Influence of k in KNN method ####

## Setting up parameters of estimation
ns=c(2000)
sds=c(1)
i_sd=sds

if (analysis_type=="short"){
  knns=c(2,20,50)
} else {
  knns=c(2,4,6,8,10,15,20,25,30,50,100)
}

corr_init=0
i_outputs=paste(c(2,10,30),"d",sep="")

i_class="2"
i_type="multigauss_knns"

path_output_main=paste('plotsSI/synthetic_comparison/',i_type,sep="")
dir.create(path_output_main,recursive = TRUE)

OutputList=list()
OutputListT=list()

# 2.8 Setting up the specifics of the channel
classes_num = as.numeric(i_class) 
i_output="1d"
dimension_num=as.numeric(str_match(i_output,"[0-9]+"))

mean_diff_init=c(2,rep(0,(dimension_num-1) ))
mean_init=rep(0,dimension_num)
sd_init=rep(1,dimension_num)

example_means1 = mean_init
example_means2 = mean_init+mean_diff_init

example_sigma1 = (i_sd)*diag(sd_init)
example_sigma2 = (i_sd)*diag(sd_init)

#Exact solution
tempoutput_exact <- cc_exact(type="gauss",classes=i_class,
                             dimension= dimension_num,
                             par=list(mu1=example_means1,
                                      mu2=example_means2,
                                      sd1=example_sigma1,
                                      sd2=example_sigma2),
                             x_length=10000,s_length=200)

for (i_output in i_outputs) {
  
  classes_num = as.numeric(i_class) 
  dimension_num=as.numeric(str_match(i_output,"[0-9]+"))
  
  mean_diff_init=c(2,rep(0,(dimension_num-1) ))
  mean_init=rep(0,dimension_num)
  sd_init=rep(1,dimension_num)
  
  example_means1 = mean_init
  example_means2 = mean_init+mean_diff_init
  
  example_sigma1 = (i_sd)^2*diag(sd_init)
  example_sigma2 = (i_sd)^2*diag(sd_init)
  
  OutputList[[as.character(i_output)]]=list()
  OutputListT[[as.character(i_output)]]=list()
  
  
  for (i_n in ns){
    print(paste(i_n,i_output,sep="_"))
    path_output_single=paste(path_output_main,
                             'dim_',dimension_num,'_N_',i_n,'/',
                             sep="")
    
    OutputList[[as.character(i_output)]][[as.character(i_n)]]=list()  
    OutputListT[[as.character(i_output)]][[as.character(i_n)]]=list()
    
    OutputList[[as.character(i_output)]][[as.character(i_n)]]$exact=tempoutput_exact
    
    #LogReg Method
    t1=proc.time()
    cl=parallel::makeCluster(boot_cores)
    doParallel::registerDoParallel(cl)
    OutList=foreach::foreach(i=1:boot_num,.packages = c("SLEMI","mvtnorm",
                                                        "FNN","TDA","DEoptimR")) %do% {
      
      tempdata = data.frame(signal = c(t(replicate(i_n,LETTERS[1:classes_num]))) ,
                            rbind(
                              rmvnorm(i_n,mean=example_means1,sigma=example_sigma1),
                              rmvnorm(i_n,mean=example_means2,sigma=example_sigma2) 
                            ))
      
      path_output_single_boot=paste(path_output_single,i,"/",sep="")
      dir.create(path_output_single_boot,recursive = TRUE)
      tempoutput  <- capacity_logreg_main(dataRaw=tempdata,
                                          signal="signal", response=paste("X",1:dimension_num,sep=""),
                                          formula_string = paste("signal~", paste(paste("X",1:dimension_num,sep=""),collapse="+"),sep=""),
                                          output_path=path_output_single_boot,plot_height=8,plot_width=12,
                                          testing=FALSE,model_out = FALSE,data_out = FALSE)
      tempout=tempoutput$cc
      rm(tempdata,tempoutput)
      tempout
    }
    stopCluster(cl)
    OutputListT[[as.character(i_output)]][[as.character(i_n)]]$lr= (proc.time()-t1)[3]
    OutputList[[as.character(i_output)]][[as.character(i_n)]]$lr=OutList 
    
    OutputListT[[as.character(i_output)]][[as.character(i_n)]]$knn=list()
    OutputList[[as.character(i_output)]][[as.character(i_n)]]$knn=list()
    
    # KNN method
    for (i_knn in knns){
      
      t1=proc.time()
      cl=parallel::makeCluster(boot_cores)
      doParallel::registerDoParallel(cl)
      OutList=foreach::foreach(i=1:boot_num,.packages = c("SLEMI","mvtnorm",
                                                          "FNN","TDA","DEoptimR")) %do% {
        
        tempdata = data.frame(signal = c(t(replicate(i_n,LETTERS[1:classes_num]))) ,
                              rbind(
                                rmvnorm(i_n,mean=example_means1,sigma=example_sigma1),
                                rmvnorm(i_n,mean=example_means2,sigma=example_sigma2) 
                              ))
        
        ni_class=as.numeric(i_class)
        
        knndens=list()
        for (is in LETTERS[1:classes_num]){
          knndens[[is]]=knnDE(tempdata[tempdata$signal==is,-1],tempdata[,-1],i_knn)
        }
        
        eval_f_args<-function(q,data,knndens){
          
          signals=as.character(unique(data[,1]))
          fri=lapply(signals,function(x) { knndens[[x]][data[,1]==x]  } )
          fr=Reduce("+",mapply(function(x,y){ y*knndens[[x]] },signals,q,SIMPLIFY=FALSE))
          
          hrs=-sum( mapply(function(x,y){
            sum(log2(x))*y*(1/length(x))
          },fri,q,SIMPLIFY=TRUE)  )
          
          hr=-sum( mapply(function(x,y){
            sum(log2(fr[data[1]==x]))*y*(1/length(x[data[1]==x]))
          },signals,q,SIMPLIFY=TRUE)  )
          
          return(hrs-hr)
        }
        
        
        eval_f<-function(q){return(eval_f_args(q,data=tempdata,knndens= knndens))}
        eval_g<-function(q){return(1-sum(q))}
        
        resDE<-JDEoptim(rep(0,ni_class),rep(1,ni_class),
                        fn = eval_f, constr = eval_g, meq = 1,
                        tol = 1e-6, trace = FALSE, triter = 50,maxiter = 6000)
        
        tempout=-resDE$value
        
        rm(tempdata,knndens)
        tempout
      }
      stopCluster(cl)
      OutputListT[[as.character(i_output)]][[as.character(i_n)]]$knn[[as.character(i_knn)]]=(proc.time()-t1)[3]
      OutputList[[as.character(i_output)]][[as.character(i_n)]]$knn[[as.character(i_knn)]]=OutList 
      
    }
  }
}
saveRDS(OutputListT,file=paste0(path_output_main,"/times.rds")) #saving objects in files
saveRDS(OutputList,file=paste0(path_output_main,"/knn.rds"))



# 2.9 Generating graphs
df=melt(OutputListT)
df$L4=as.numeric(as.character(df$L4))
df$L1=factor(df$L1,levels=paste(c(2,10,30),"d",sep=""))
df$L2=factor(df$L2,levels=c(2000))
plot=ggplot(data=df,aes(x=L4,y=value))+geom_line()+facet_grid(L2~L1)
ggsave(plot,file="plotsSI/plot_S1D.pdf",height=6,width=8)


dff_lr=melt(lapply(OutputList,function(x) lapply(x,function(y) lapply(y[c("lr")], function(z) z  ) ) ))
dff_lr$L5=NA
dff_knn=melt(lapply(OutputList,function(x) lapply(x,function(y) lapply(y[c("knn")], function(z) lapply(z,function(w) w)  ) ) ))
dff_ex=melt(lapply(OutputList,function(x) lapply(x,function(y) lapply(y[c("exact")], function(z) z$cc  ) ) ))
dff_ex$L4=NA
dff_ex$L5=NA

dff=rbind(dff_lr,dff_knn,dff_ex)

dff$L1=factor(dff$L1,levels=paste(c(2,10,30),"d",sep=""))
dff$L2=as.numeric(as.character(dff$L2))
dff$L4=as.numeric(as.character(dff$L4))

plot=ggplot(data=dff[dff$L3=="knn",],aes(x=L4,y=value))+facet_grid(L1~L2)+stat_summary()+
  geom_hline(yintercept = mean(dff_ex$value))+
  scale_x_continuous("k in knn")+scale_y_continuous("capacity")+
  theme_publ(version=2)
  
ggsave(plot,file="plotsSI/plot_S1E.pdf",height=6,width=8)



#### 3. Validation ####
## After running this chunk, plots similar to those in Fig. S1 will be created in directory plotsSI/synthetic_validation/
## Without any changes, it will take ~1hour to re-calculate all examples

full_path_out = paste('plotsSI/synthetic_validation/',sep="")

# 3.1 Set parameters of the estimation
ns=c(500,1000,2000) # sample sizes of synthetic data

if (analysis_type=="short"){
    sds=c(0.01,1,10)  # standard deviations of the example
  } else {
    sds=c(0.01,0.1,0.25,0.5,1,2,4,5,10)
  }


mean_diff_init=1 #auxillary variables for mean specification
mean_diff_step=0.1
mean_init=0
sd_init=1



## 3.2 Estimate capacity of Exponential channel ####
print("Exponential Output: starting")
i_output="1d"
i_class="5"
i_type="exp"

path_output_main=paste(full_path_out,i_type,'/',sep="")
dir.create(path_output_main,recursive = TRUE)

OutputList=list()
classes_num = as.numeric(i_class) 

for (i_var in sds){
  print(paste(i_type,i_var,sep="_"))
  example_means =  seq(from= i_var, by= 1,length=classes_num)
  tempoutput_exact <- cc_exact(type="exp",classes=i_class,dimension=str_match(i_output,"[0-9]+"),
                               par=list(lam=1/example_means),
                               x_length=500,mc_length=10000,mc_cores=boot_cores) #calculate true channel capacity
  OutputList[[as.character(i_var)]]=list()
  OutputList[[as.character(i_var)]]$exact=tempoutput_exact
  for (i_n in ns){
    path_output_single=paste(path_output_main,
                             'sd_',i_var,'_N_',i_n,'/',
                             sep="")
    
    cl=parallel::makeCluster(boot_cores)
    doParallel::registerDoParallel(cl)
    OutList=foreach::foreach(i=1:boot_num,.packages = c("SLEMI")) %dopar% {
      tempdata = data.frame(signal = c(t(replicate(i_n,LETTERS[1:classes_num]))) ,
                            output = c(matrix(rexp(classes_num*i_n,rate=1/example_means),ncol=classes_num,byrow=TRUE) ))
      path_output_single_boot=paste(path_output_single,i,"/",sep="")
      dir.create(path_output_single_boot,recursive = TRUE)
      tempoutput  <- capacity_logreg_main(dataRaw=tempdata,
                                          signal="signal", response="output",
                                          formula_string = "signal~output",
                                          output_path=path_output_single,
                                          testing=FALSE,model_out=FALSE,data_out=FALSE)
      tempoutput
    }
    parallel::stopCluster(cl)
    OutputList[[as.character(i_var)]][[as.character(i_n)]]=OutList
  }
}
print("Exponential Output: finished")
saveRDS(OutputList,file=paste(path_output_main,"outputdata.rds",sep="")) # saving in a file





## 3.3 Estimate capacity of Gamma channel####
print("Gamma Output: starting")
i_output="1d"
i_class="5"
i_type="gamma"

path_output_main=paste(full_path_out,i_type,'/',sep="")
dir.create(path_output_main,recursive = TRUE)

OutputList=list()
classes_num = as.numeric(i_class) 

example_means = mean_init + cumsum( seq(from= mean_diff_init, to= mean_diff_init,length=classes_num) ) 

for (i_var in sds){
  print(paste(i_type,i_var,sep="_"))
  example_shape =  (example_means)^2/(i_var^2)
  example_scale = (i_var^2)/example_means
  
  tempoutput_exact <- cc_exact(type="gamma",classes=i_class,dimension=str_match(i_output,"[0-9]+"),
                               par=list(shape=example_shape,scale=example_scale),
                               x_length=1000,mc_length=10000,mc_cores=boot_cores) #calculate true channel capacity
  OutputList[[as.character(i_var)]]=list()
  OutputList[[as.character(i_var)]]$exact=tempoutput_exact
  for (i_n in ns){
    
    path_output_single=paste(path_output_main,
                             'sd_',i_var,'_N_',i_n,'/',
                             sep="")
    
    cl=parallel::makeCluster(boot_cores)
    doParallel::registerDoParallel(cl)
    OutList=foreach::foreach(i=1:boot_num,.packages = c("SLEMI"),.export=c("x_log_y")) %dopar% {
      tempdata = data.frame(signal = c(t(replicate(i_n,LETTERS[1:classes_num]))) ,
                            output = c(matrix(rgamma(classes_num*i_n,shape=example_shape,scale=example_scale),ncol=classes_num,byrow=TRUE) ))
      path_output_single_boot=paste(path_output_single,i,"/",sep="")
      dir.create(path_output_single_boot,recursive = TRUE)
      tempoutput  <- capacity_logreg_main(dataRaw=tempdata,
                                          signal="signal", response="output",
                                          formula_string = "signal~output+I(output^2)+I(log(output+1e-16))",
                                          output_path=path_output_single,
                                          testing=FALSE,model_out=FALSE,data_out=FALSE,scale=FALSE)
      tempoutput
    }
    parallel::stopCluster(cl)
    OutputList[[as.character(i_var)]][[as.character(i_n)]]=OutList
  }
}
saveRDS(OutputList,file=paste(path_output_main,"outputdata.rds",sep=""))  # saving in a file
print("Gamma Output: finished")




## 3.4 Estimate capacity of Lognormal channel####
print("LogNormal Output: starting")
i_output="1d"
i_class="5"
i_type="lnorm"

path_output_main=paste(full_path_out,i_type,'/',sep="")
dir.create(path_output_main,recursive = TRUE)

OutputList=list()
classes_num = as.numeric(i_class) 

example_means = (mean_init + cumsum( seq(from= mean_diff_init, to= mean_diff_init,length=classes_num) ) )

for (i_var in sds){
  print(paste(i_type,i_var,sep="_"))
  example_sds = rep( sd_init*i_var ,classes_num)
  example_logmeans=log(example_means/sqrt(1+ (example_sds^2/example_means^2) ))
  example_logsds=sqrt(log(1+(example_sds^2/example_means^2)))
  
  tempoutput_exact <- cc_exact(type="lnorm",classes=i_class,dimension=str_match(i_output,"[0-9]+"),
                               par=list(mu=example_logmeans,sd=example_logsds),
                               x_length=2000,mc_length=10000,mc_cores=boot_cores) #calculate true channel capacity
  OutputList[[as.character(i_var)]]=list()
  OutputList[[as.character(i_var)]]$exact=tempoutput_exact
  for (i_n in ns){
    
    path_output_single=paste(path_output_main,
                             'sd_',i_var,'_N_',i_n,'/',
                             sep="")
    
    cl=parallel::makeCluster(boot_cores)
    doParallel::registerDoParallel(cl)
    OutList=foreach::foreach(i=1:boot_num,.packages = c("SLEMI"),.export=c("x_log_y")) %dopar% {
      tempdata = data.frame(signal = c(t(replicate(i_n,LETTERS[1:classes_num]))) ,
                            output = c(matrix(rlnorm(classes_num*i_n,meanlog =example_logmeans ,sdlog = example_logsds),ncol=classes_num,byrow=TRUE) ))
      path_output_single_boot=paste(path_output_single,i,"/",sep="")
      dir.create(path_output_single_boot,recursive = TRUE)
      tempoutput  <- capacity_logreg_main(dataRaw=tempdata,
                                          signal="signal", response="output",
                                          formula_string = "signal~output+I(log(output+1))+I((log(output+1))^2 )+I((log(output+1))^3)+I(sqrt(output))",
                                          output_path=path_output_single,
                                          testing=FALSE,model_out=FALSE,data_out=FALSE,scale=FALSE)
      tempoutput
    }
    parallel::stopCluster(cl)
    OutputList[[as.character(i_var)]][[as.character(i_n)]]=OutList
    
  }
}
saveRDS(OutputList,file=paste(path_output_main,"outputdata.rds",sep="")) #saving results in a file
print("LogNormal Output: finished")




## 3.5 Estimate capacity of Gaussian channel ####
print("Normal Output: starting")
i_output="1d"
i_class="5"
i_type="norm"
path_output_main=paste(full_path_out,i_type,'/',sep="")
dir.create(path_output_main,recursive = TRUE)

OutputList=list()
classes_num = as.numeric(i_class) 

example_means = mean_init + cumsum( seq(from= mean_diff_init, to= mean_diff_init,length=classes_num) ) 

for (i_var in sds){
  print(paste(i_type,i_var,sep="_"))
  example_sds = rep( sd_init*i_var ,classes_num)
  tempoutput_exact <- cc_exact(type="gauss",classes=i_class,dimension=str_match(i_output,"[0-9]+"),
                               par=list(mu=example_means,sd=example_sds),
                               x_length=500,mc_length=10000,mc_cores=boot_cores) #calculate true channel capacity
  OutputList[[as.character(i_var)]]=list()
  OutputList[[as.character(i_var)]]$exact=tempoutput_exact
  for (i_n in ns){
    
    path_output_single=paste(path_output_main,
                             'sd_',i_var,'_N_',i_n,'/',
                             sep="")
    
    cl=parallel::makeCluster(boot_cores)
    doParallel::registerDoParallel(cl)
    OutList=foreach::foreach(i=1:boot_num,.packages = c("SLEMI"),.export=c("x_log_y")) %dopar% {
      tempdata = data.frame(signal = c(t(replicate(i_n,LETTERS[1:classes_num]))) ,
                            output = c(matrix(rnorm(classes_num*i_n,mean=example_means,sd=example_sds),ncol=classes_num,byrow=TRUE) ))
      path_output_single_boot=paste(path_output_single,i,"/",sep="")
      dir.create(path_output_single_boot,recursive = TRUE)
      tempoutput  <- capacity_logreg_main(dataRaw=tempdata,
                                          signal="signal", response="output",
                                          formula_string = "signal~output",
                                          output_path=path_output_single,
                                          testing=FALSE,model_out=FALSE,data_out=FALSE,scale=FALSE)
      tempoutput
    }
    parallel::stopCluster(cl)
    OutputList[[as.character(i_var)]][[as.character(i_n)]]=OutList
    
  }
}
saveRDS(OutputList,file=paste(path_output_main,"outputdata.rds",sep="")) #saving results in a file
print("Normal Output: finished")



## 3.6 Generating plots
dataE=readRDS(file=paste0(full_path_out,"exp/outputdata.rds"))
dataG=readRDS(file=paste0(full_path_out,"gamma/outputdata.rds"))
dataN=readRDS(file=paste0(full_path_out,"norm/outputdata.rds"))
dataLN=readRDS(file=paste0(full_path_out,"lnorm/outputdata.rds"))
dataN$`0.01`$exact$cc=log2(5)
dataLN$`0.01`$exact$cc=log2(5)

temppath=paste("plotsSI/synthetic_validation/plots_exp/",sep="")
dir.create(temppath,recursive = TRUE)
plot12=tempFun3_plot(dataE,5,"1d",temppath,plot_height=6,plot_width=10,ylimits=c(0,NA),sdlimit=10)

temppath=paste("plotsSI/synthetic_validation/plots_gamma/",sep="")
dir.create(temppath,recursive = TRUE)
plot22=tempFun3_plot(dataG,5,"1d",temppath,plot_height=6,plot_width=10,ylimits=c(0,NA),sdlimit=10)

temppath=paste("plotsSI/synthetic_validation/plots_normal/",sep="")
dir.create(temppath,recursive = TRUE)
plot32=tempFun3_plot(dataN,5,"1d",temppath,plot_height=6,plot_width=10,ylimits=c(0,NA),sdlimit=10)

temppath=paste("plotsSI/synthetic_validation/plots_lognormal/",sep="")
dir.create(temppath,recursive = TRUE)
plot42=tempFun3_plot(dataLN,5,"1d",temppath,plot_height=6,plot_width=10,ylimits=c(0,NA),sdlimit=10)

plot=grid.arrange(plot12,plot22,plot32,plot42,ncol=2,nrow=2)
ggsave(plot,file=paste("plotsSI/plot_S2.pdf",sep=""),height=12,width=20)


