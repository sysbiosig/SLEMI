#' Examplary data set I
#'
#' A dataset describing simple one dimensional intput - one dimensional output channel 
#' with 500 observations per input. Each observation can be in addition of 3 types
#' that occurs with propensities (0.6,0.3,0.1), resepctively. 
#' Conditional output distributions are Gaussian.
#' 
#' @format A data frame with 1500 rows and 3 variables:
#' \describe{
#' \item{signal}{Label of input}
#' \item{response}{The value of output}
#' \item{sideVar}{Label of the type of given observation}
#' }
#' @source synthetic
"data_example1"

#' Examplary data set II
#'
#' A dataset describing a channel with 3 possible inputs and 3-dimensional output
#' with 500 observations per input.
#' Conditional output distributions are multivariate Gaussians.
#' 
#' @format A data frame with 1500 rows and 4 variables:
#' \describe{
#' \item{signal}{Label of input}
#' \item{X1}{The value of first dimension of output}
#' \item{X2}{The value of second dimension of output}
#' \item{X3}{The value of third dimension of output}
#' }
#' @source synthetic
"data_example2"


#' Data from experiment with NFkB pathway
#'
#' In the paper describing methodological aspects of our algorithm we present the analysis of information transmission 
#' in NfkB pathway upn the stimulation of TNF-$\alpha$. Experimental data from this experiment in the form of single-cell 
#' time series are attached to the package as a data.frame object and can be accessed using `data_nfkb` variable.

#' Each row of `data_nfkb` represents a single observation of a cell. Column 'signal' indicates the level of TNF-$\alpha$ 
#' stimulation for a given cell, while columns 'response_T', gives the normalised ratio of nuclear and cytoplasmic transcription 
#' factor as described in Supplementary Methods of the corresponding publication. 
#' 
#' For each concentration, there are at least 1000 single-cell observation (with the exception of 0.5ng stimulation, 
#' where the number of identified cells is almost 900)
#' @format A data frame with 15632 rows and 42 variables:
#' \describe{
#' \item{signal}{Level of TNFa stimulation}
#' \item{response_T}{The concentration of normalised NfkB transcription factor, measured at time T}
#' #' }
#' @source in-house experimental data
"data_nfkb"

