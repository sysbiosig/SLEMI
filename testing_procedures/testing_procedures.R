### ### ### ### ### ### ### ### ### ### ### ### 
###                                         ###
### title: Script for testing procedures    ###
###         of SLEMI R package              ###
### author: T. Jetka, K. Nienaltowski,      ###
###         T. Winarski, M. Komorowski      ###
### date: 22.07.2018                        ###
### Submission of manuscript:               ###
###                                         ###
### Information-theoretic analysis of       ###
### multivariate signaling responses        ###
### using SLEMI                             ###
###                                         ###
### ### ### ### ### ### ### ### ### ### ### ### 

#### Descritption ####
### ### ###
### Script, that demonstrates how to use and test SLEMI package 
### to calculate channel capacity from experimental data. 
### A full analysis starting from 
### - preparation of data 
### - estimation of channel capacity 
### - visualisation of results and diagnostic statistics 
###
### Script also verify other functionalities accompanying
### the package will be presented.
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
require("mvtnorm")   # install.packages("mvtnorm")
require("ggplot2")   # install.packages("ggplot2")
require("gridExtra") # install.packages("gridExtra")

# FOR USER: install tested package SLEMI
#install_github("sysbiosig/SLEMI")


#### Script configuration ####
### 1. 
# FOR USER: set working directory to folder with script
path_wd <- "~/path/to/folder/testing_procedures/" 
dir.create(path_wd,recursive = TRUE)
setwd(path_wd)

### 2. 
# Load SLEMI package 
library(SLEMI)

### 3. 
# Set a seed of a random number generator for reproducability 
set.seed(3349)

### 4. 
# FOR USER: set number of cores
cores_num <- 1

### 5. 
# FOR USER: variable that determines whether to display graphs
display_plots <- TRUE

#### Preparing data ####
tempdata <- data_nfkb[,1:4]
tempdata <- tempdata[!apply(tempdata,1,function(x) any(is.na(x))),]
tempdata <- rbind(tempdata[1:4,],
                  tempdata[10001:10004,],
                  tail(tempdata,4))
tempdata <- tempdata[order(tempdata$signal),]
print(head(tempdata))
#       signal response_0 response_3 response_6
# 1        0ng  0.3840744  0.4252835  0.4271986
# 2        0ng  0.4709216  0.5777821  0.5361948
# 3        0ng  0.4274474  0.6696011  0.8544916
# 4        0ng  0.4776072  0.4409360  0.4224470
# 16229  100ng  0.6236179  0.6154709  1.0816505
# 16240  100ng  1.3534083  3.0158004  5.1592848

#### Running simple examples : Generating synthetic dataset ####
# Channel with a Gaussian output 

### 6. Set parameters of the example
# sample size for each input concentration; 
# change to investigate its influence on estimation
n_sample <- c(1000) 
# standard deviation of Y|X=x;
# change to investigate its influence on estimation
dist_sd  <- 1 
# number of concentration of input X considered;
# change to investigate its influence on estimation
input_num <- 6 

### 7. Create output directory
i_type <- "testing_basic" 
path_output_main <- 
  paste('output/',
        i_type,'/',
        sep="")
dir.create(
  path_output_main,
  recursive = TRUE)

### 8. Generating synthetic data
# concentration of input; spans from 0 to saturation.
xx <- 
  signif(
    c(0,
      exp(
        seq(
          from = log(0.01),
          to = log(100),
          length.out = input_num-1))),
    digits = 2) 

# mean of the dose-response relation; Michealis Menten assumed
example_means <- 10*(xx/(1+xx)) 
example_sds <- rep(dist_sd, input_num)
tempdata <-
  data.frame(
    signal = c(t(replicate(n_sample, xx))),
    output = c(matrix(
      rnorm(n = input_num*n_sample,
            mean = example_means,
            sd = example_sds),
      ncol = input_num,
      byrow = TRUE)))
tempdata$signal <-
  factor(
    x = tempdata$signal,
    levels = sort(unique(tempdata$signal)))
print(head(tempdata))
# signal            output
#      1      0  2.5174406
#      2      0  1.5522436
#      3      0  0.8514287
#      4      0 -2.3999916
#      5      0 -0.4190470
#      6      0  1.7573737

### 9 Preview the distribution of the data
g_plot <- ggplot(
  data = tempdata,
  aes(x = factor(signal),
      group = signal,
      y = output)) + 
  geom_boxplot() +
  theme_publ(version = 2)

if(display_plots){
  print(g_plot)
}

### 10. Set required parameters for the algorithm
signal_name <- "signal"
response_name <- "output"

### 11. Estimate channel capacity (takes several seconds)
tempoutput  <- 
  capacity_logreg_main(
    dataRaw = tempdata, 
    signal = signal_name,
    response = response_name, 
    output_path = path_output_main
  )

### 12 Results of the estimation
# channel capacity in bits
print(paste("Channel Capacity (bit):", 
            tempoutput$cc,
            sep=" ")) 
# [1] "Channel Capacity (bit): 1.59504052094311"

# optimal probabilities
print(paste("Optimal input probabilities:",
            paste(tempoutput$p_opt,
                  collapse = ", "),
            sep = " " )) 
# [1] "Optimal input probabilities: 0.216995711379725, 0.00128924561759806, 
# 0.13292998209147, 0.30468664259557, 0.14246143408601, 0.201636984229627"

# accuracy of classification
print(paste("Accuracy of classification:",
            tempoutput$regression$overall[1],
            sep = " " )) 
# [1] "Accuracy of classification: 0.601"

# time of computations
print(paste("Time of computations (sec.):",
            tempoutput$time[3],
            sep = " " )) 
# [1] "Time of computations (sec.): 5.68300000000002"

### 13 Inspect object generated during computation
# verification if results saved in rds are the same;
# full output of the estimation should be saved in  
# output_path/output.rds;
tempoutput_rds = readRDS(
  paste0(
    path_output_main,
    "/output.rds"))

# assert channel capacity in bits
print(paste("Channel Capacity assertion:",
            tempoutput_rds$cc == tempoutput$cc))
# [1] "Channel Capacity assertion: TRUE"

# assert optimal probabilities
print(paste("Optimal input probabilities assertion:", 
            sum(tempoutput_rds$p_opt != tempoutput$p_opt) == 0)) 
# [1] "Optimal input probabilities assertion: TRUE"

# assert accuracy of classification
print(paste("Accuracy of classification assertion:",
            tempoutput_rds$regression$overall[1] == 
              tempoutput$regression$overall[1]))
# [1] "Accuracy of classification assertion: TRUE"

# assert time of computations
print(paste("Time of computations assertion:",
            tempoutput_rds$time[3] == 
              tempoutput$time[3]))
# [1] "Time of computations assertion: TRUE"


### 14 Explore visualisation of the results. 
# See pdf files in output_path/ for the visualisation of the data and 
# capacity estimation. 
# - "MainPlot.pdf" - presents the most important information: mean,
# input-output relation, distributions of output and channel capacity (see below).
# - graphs, as gg or gtable objects, are saved in `logGraphs` element of function output.

if(display_plots){
  grid.arrange(tempoutput$logGraphs[[9]])
}

### 15 Generate graphs of different size 
# specify paramters `plot_width` and `plot_height`
plot_width_new  <- 10
plot_height_new <- 8 
i_type <- "testing_basic_graphs_size"
path_output_main <- paste('output/', i_type,'/', sep="")
dir.create(path_output_main,
           recursive = TRUE)
tempoutput  <- capacity_logreg_main(
  dataRaw = tempdata, 
  signal = signal_name, 
  response = response_name, 
  output_path = path_output_main,
  plot_width = plot_width_new,
  plot_height = plot_height_new 
)

### 16 Run analysis without generating graphs
# set argument `graphs` to FALSE
graphs_generate <- FALSE
i_type <- "testing_basic_nographs"
path_output_main <- paste('output/',i_type,'/',sep="")
dir.create(path_output_main,
           recursive = TRUE)
tempoutput  <- capacity_logreg_main(
  dataRaw = tempdata,
  signal = signal_name, 
  response = response_name, 
  output_path = path_output_main,
  graphs = graphs_generate
)

### 17 Run analysis with the minimal output 
# set `graphs`, `scale`, `dataout` and `model_out` to FALSE; 
# respectively, this prevents creating graphs, scaling of the
# data, including data and regression model in returned list.
graphs_generate <- FALSE
data_rescale <- FALSE
data_save <- FALSE
model_save <- FALSE
i_type <- "testing_basic_minimalOutput"
path_output_main=paste('output/',i_type,'/',sep="")
dir.create(path_output_main,recursive = TRUE)
tempoutput  <- capacity_logreg_main(
  dataRaw = tempdata,
  signal = signal_name, 
  response = response_name, 
  output_path = path_output_main,
  graphs = graphs_generate,
  scale = data_rescale,
  dataout = data_save,
  model_out = model_save
)

#### Diagnostic tests ####

### 18  Set parameters of diagnostic tests
# Set a seed of a random number generator for reproducability 
seed_to_use <- 12345
# number of bootstrap repetitions 
bootstrap_num <- 20
# fraction of data to sample
bootstrap_frac <- 0.8 
# number of repetition of overfitting test
overfitting_num <- 20 
# fraction of data to use as training sample
training_frac <- 0.6 
 
i_type <- "testing_basic_diagnostic"
path_output_main <- paste('output/',i_type,'/',sep="")
dir.create(path_output_main,
           recursive = TRUE)

### 19 Run estimation with full diagnostics 
# set `testing` argument to TRUE (takes ~XX min)

tempoutput  <- capacity_logreg_main(
  dataRaw = tempdata, 
  signal = signal_name,
  response = response_name, 
  output_path = path_output_main,
  testing = TRUE,
  plot_width = 10 ,
  plot_height = 8,
  TestingSeed = seed_to_use,
  testing_cores = cores_num,
  boot_num = bootstrap_num,
  boot_prob = bootstrap_frac,
  traintest_num = overfitting_num,
  partition_trainfrac = training_frac
)

### 20 Inspect the results of diagnostic tests. 
# Results of diagnostic tests is saved in `testing` element of output list

#  channel capacity in bits
print(paste("Channel Capacity, bootstrap mean (sd): ",
            mean(sapply(tempoutput$testing$bootstrap,
                        function(x) x$cc)),
            "(",sd(sapply(tempoutput$testing$bootstrap,
                          function(x) x$cc)),")",
            sep = "" )) 
# [1] "Channel Capacity, bootstrap mean (sd): 1.59379723093333(0.00408538627026014)"

#  time of computations
print(paste("Time of computations (sec.):",
            tempoutput$time[3],
            sep = " ")) 
# [1] "Time of computations (sec.): 32.939"

### 21 Visualistion of diagnostic tests 
# See the visualisation in the output directory - MainPlot.pdf.
# For each diagnostic test, there is a corresponding histogram of
# calculated capacities. 
# This graph is also obtainable from the functions results.

if(display_plots){
  grid.arrange(tempoutput$logGraphs[[9]])
}

### 22 p-values of diagnostic tests 
# For each diagnostic test, we provide left- and right-tailed empirical 
# p-values of obtained channel capacity. They indicated if the regime 
# of data bootstraping or dividing data into training and testing sample 
# influence the calculation of capacity in a significany way. 
# A small p-value in any of these tests (e.g. <0.05) means a problem with
# stability of channel capacity estimation and a possible bias due to too
# small sample size. P-values are printed on the MainPlot.pdf graph
# or can be obtained in
print(paste("P-values:",
            tempoutput$testing_pv$bootstrap[1], 
            " (left-tailed); ",
            tempoutput$testing_pv$bootstrap[2], 
            " (right-tailed) ",
            sep = ""))
# [1] "P-values:0.7 (left-tailed); 0.3 (right-tailed) "

### 23 summary of diagnostic tests
print(summary(tempoutput))
#            Length Class           Mode   
# regression  6     confusionMatrix list   
# model      26     multinom        list   
# p_opt       6     -none-          numeric
# cc          1     -none-          numeric
# testing     2     -none-          list   
# testing_pv  2     -none-          list   
# time        5     proc_time       numeric
# params      3     -none-          numeric
# data        2     data.frame      list   
# logGraphs   9     -none-          list   
#### Channel with a multi-dimensional output ####

### 24. Load synthetic dataset. 
# We created a synthetic data set, where output is three dimensional
# and there are three different inputs.
# It accompanies our package in variable data example2
tempdata <- data_example2

# inspect data
head(tempdata)
#   signal          X1          X2          X3
# 1      0  0.15580917  0.18045768  0.02999694
# 2      0 -0.16423212  0.04744891 -0.42674054
# 3      0  0.13896586 -0.06779270 -0.05983137
# 4      0 -0.16338335  0.04204709  0.40134308
# 5      0  0.07555051  0.09766543 -0.14256936
# 6      0  0.16996827 -0.21497172 -0.09494647

# 25. Channel capacity of synthetic data
# preapration of channel capacity computation
signal_name <- "signal"
response_name <- c("X1","X2","X3")
i_type <- "testing_multivariate" 
path_output_main <- paste('output/',
                          i_type,
                          '/',
                          sep = "")
dir.create(path_output_main,
           recursive = TRUE)

# computation of channel capacity
tempoutput  <- capacity_logreg_main(
  dataRaw = tempdata, 
  signal = signal_name,
  response = response_name, 
  output_path = path_output_main
)

### 26 Print output of the estimation in the console
# channel capacity in bits
print(paste("Channel Capacity (bit):",
            tempoutput$cc,
            sep = " ")) 
# [1] "Channel Capacity (bit): 1.58496240473138"

# optimal probabilities
print(paste("Optimal input probabilities:",
            paste( 
              tempoutput$p_opt,
              collapse = ", "),
            sep = " ")) 
# [1] "Optimal input probabilities: 0.333333327033102, 0.333333320455409,
# 0.333333352511489"

# accuracy of classification
print(paste("Accuracy of classification:",
            tempoutput$regression$overall[1],
            sep = " ")) 
# [1] "Accuracy of classification: 1"

# time of computations
print(paste("Time of computations (sec.):",
            tempoutput$time[3],
            sep = " ")) 
# [1] "Time of computations (sec.): 1.25699999999983"

### 27. Visualistion of analysis
# See the visualisation of results in  MainPlot.pdf 
# in output directory or in object

if(display_plots){
  grid.arrange(tempoutput$logGraphs[[9]])
}
