# Statistical Learning based Estimation of Mutual Information (SLEMI)
 The R package SLEMI is designed to estimate channel capacity between finite state input and multidimensional output from experimental data. For efficient computations, it uses iterative algorithm based on logistic regression.  The core function `capacity_logreg_main()` is the basic interface to all functionalities provided in the package. A comprehensive documentation is available in directory [`vignette/SLEMI_vignette.pdf`](https://github.com/sysbiosig/SLEMI/blob/master/vignette/SLEMI_vignette.pdf).
 ## Setup
 ### Requirements - Hardware
  + A 32 or 64 bit processor (recommended: 64bit)
  + 1GHz processor (recommended: multicore for a comprehensive analysis)
  + 2GB MB RAM (recommended: 4GB+, depends on the size of experimental data)
 ### Requirements - Software
The main software requirement is the installation of the R environment (version: >= 3.6), which can be downloaded from [R project website](https://www.r-project.org) and is distributed for all common operating systems.  We tested the package in R environment installed on Windows 7, 10; Mac OS X 10.11 - 10.13 and Ubuntu 18.04 with no significant differences in the performance. The use of a dedicated Integrated development environment (IDE), e.g. [RStudio](https://www.rstudio.com) is recommended. 
 Apart from base installation of R, SLEMI requires following packages:
 1. for installation 
   + devtools
  
2. for estimation
  
  + e1071
  + Hmisc
  + nnet
  + glmnet
  + caret
  + doParallel (if parallel computation are needed)
  
3. for visualisation
  + ggplot2
  + ggthemes
  + gridExtra
  + corrplot
  
4. for data handling
  + reshape2
  + stringr
  + plyr
 Each of the above packages can be installed by running 
> `install.packages("name_of_a_package")`
 in the R console.
 ### Installation
 In order to install the package use following commands in R's console
 > `# install.packages("devtools") # run if not installed`
 
 > `library(devtools)`
 
 > `install_github("sysbiosig/SLEMI")`
 
 ## Basic usage
 The package is based on a main wrapper function - `capacity_logreg_main()` for calculation of channel capacity, which calls specific methods implemented within this package. Similarly, functions `mi_logreg_main()` can be used to estimate mutual information, while `prob_discr_pairwise()` to compute probabilities of discrimination between two different input states.
 ### Preaparing data
 For the calculation of channel capacity between X and Y you need structure experimental data into a single `data.frame` object with observations in rows, one column with values of input (X), prefferably of `factor` type and columns with measured output (Y) of `numeric` type.
 ### Run
 In order to estimate channel capacity, using basic logistic regression model, call
> `capacity_logreg_main(dataRaw, signal, response, output_path)`
where: 
* `dataRaw` is a data.frame with experimental data as described above
* `signal` is a character indicating the name of column in `dataRaw` with the input (X)
* `response` is a character vector indicating names of columns in `dataRaw` with output (Y) variables
* `output_path` is a character with the directory, to which results of the estimation should be saved
 ### Results
 The function `capacity_logreg_main` returns a list, whose main elements are
 * cc - channel capacity estimate (in bits)
* p_opt - numeric bector with the optimal input distribution
* model - `nnet` object describing fitted logistic regression model
 For convenience of further analysis, this list is saved in `output_path` directory in a file `output.rds`. In addition to that, a set of exploratory graphs are created to visualise obtained estimates.
 ## Examples
 Additional examples of using package with some background on information theory is given in [`paper/TestingProcedures.pdf`](https://github.com/sysbiosig/SLEMI/blob/master/paper/TestingProcedures.pdf) and implemented in script [`paper/testing_procedures.R`](https://github.com/sysbiosig/SLEMI/blob/master/paper/testing_procedures.R). Codes used in publication are accessible from [`paper/paper_MP.R`](https://github.com/sysbiosig/SLEMI/blob/master/paper/paper_MP.R) and [`paper/paper_SI.R`](https://github.com/sysbiosig/SLEMI/blob/master/paper/paper_SI.R) respectively.
 ### Datasets
 In the manuscript describing methodological aspects of our algorithm we present the analysis of information transmission in NfkB pathway upn the stimulation of TNF-$\alpha$. Experimental data from this experiment in the form of single-cell time series are attached to the package as a data.frame object and can be accessed using `data_nfkb` variable.
 Each row of `data_nfkb` represents a single observation of a cell. Column 'signal' indicates the level of TNF-$\alpha$ stimulation for a given cell, while columns 'response_T', gives the normalised ratio of nuclear and cytoplasmic transcription factor as described in Supplementary Methods of the corresponding publication. 
 ## Other functionalities
 ### Additional paramters
 Apart from required arguments, the function `capacity_logreg_main` has also other parameters than can be used to tune the activity of the algorithm. These are
 * `model_out` (`default=TRUE`) - logical, specify if `nnet` model object should be saved into output file
* `graphs` (`default=TRUE`) - logical, controls creating diagnostic plots in the output directory.
* `plot_width` (`default = 6`) - numeric, the basic width of created plots 
* `plot_height` (`default = 4`) - numeric, the basic height of created plots
* `scale` (`default = TRUE`) - logical, value indicating if the columns of `dataRaw` are to be centered and scaled, what is usually recommended for the purpose of stability of numerical computations. From a purely theoretical perspective, such transformation does not influence the value of channel capacity.
* `lr_maxit` (`default = 1000`) - (argumnet of `nnet` package) a maximum number of iterations of optimisation step in logistic regression algorithm. Set to higher value if your data is more complex or of high dimension.
* `MaxNWts` (`default = 5000`) - (argumnet of `nnet` package) a maximum number of paramters in logistic regression model. Set to higher value if you data has many dimensions or input has many states.
 ### Diagnostic procedures
 We implemented two diagnostic procedures to control the performance of channel capacity estimation and to measure uncertainity due to finite sample size and model over-fitting. These include:
 1. Bootstrap test - capacity is re-calculated using $x$% of data, sampled from original dataset without replacement. After repeating procedure $n$ times, its standard deviation can be treated as an error of original estimate.
2. Over-fitting test - original data is divided into Training and Testing datasets. Then, logistic regression is estimated using $x$% of data (training dataset) and integrals of channel capacity are calculated via Monte Carlo using remaining $(1-x)$% of data (testing dataset). It is repeated $n$ times.
 In order to use those procedures, user must provide additional arguments to function `logreg_capacity_main()`, i.e.
 * testing (default=FALSE) - a logical value that turn on/off testing mode,
* TestingSeed (default= 1234) - the seed for the random number generator used to sample original dataset,
* testing_cores (default= 4) - a number of cores to use (via `doParallel` package) in parallel computing, 
* boot_num (default= 40) - a number of repetitions of the bootstrap , 
* boot_prob (default= 0.8) - a fraction of initial observations to  use in the bootstrap, 
* traintest_num (default= 40) - a number of repetitions of the overfitting test,
* partition_trainfrac (default= 0.6) - a fraction of initial observations to use as a training dataset in the overfitting test
 ## Support
 Please mail t.jetka at gmail.com in case of any bugs, problems and questions regarding package or inquiries regarding information theory.
 ## Reference
 Please cite
> Jetka T, Nienałtowski K, Winarski T, Błoński S, Komorowski M (2019) Information-theoretic analysis of multivariate single-cell signaling responses. PLOS Computational Biology 15(7): e1007132. https://doi.org/10.1371/journal.pcbi.1007132
 ## Licence
 SLEMI is released under the GNU licence and is freely available. A comprehensive documentation is available in directory [`vignette/SLEMI_vignette.pdf`](https://github.com/sysbiosig/SLEMI/blob/master/vignette/SLEMI_vignette.pdf).
