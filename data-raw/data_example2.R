set.seed(12345)
library(mvtnorm)

mean_diff_init=c(1,2,1.5)
mean_init=c(0,0,0)
sd_init=c(0.5,0.5,0.5)
corr_init=c(0.2,0.3,0.4)
i_n=500
i_var=0.5
classes_num = 3

example_means1 = mean_init
example_means2 = mean_init+mean_diff_init
example_means3 = mean_init+2*mean_diff_init
example_sd1=sd_init[1]*i_var
example_sd2=sd_init[2]*i_var
example_sd3=sd_init[3]*i_var
  
  t_corr_init=corr_init
  
  example_sigma1 = matrix(c(example_sd1^2,t_corr_init[1]*example_sd1*example_sd2,t_corr_init[2]*example_sd1*example_sd3,
                            t_corr_init[1]*example_sd1*example_sd2,example_sd2^2,t_corr_init[3]*example_sd2*example_sd3,
                            t_corr_init[2]*example_sd1*example_sd3,t_corr_init[3]*example_sd2*example_sd3,example_sd3^2),
                          3,3)

 # t_corr_init=corr_init
  
  example_sigma2 = matrix(c(example_sd1^2,t_corr_init[1]*example_sd1*example_sd2,t_corr_init[2]*example_sd1*example_sd3,
                            t_corr_init[1]*example_sd1*example_sd2,example_sd2^2,t_corr_init[3]*example_sd2*example_sd3,
                            t_corr_init[2]*example_sd1*example_sd3,t_corr_init[3]*example_sd2*example_sd3,example_sd3^2),
                          3,3)

 # t_corr_init=corr_init
  
  example_sigma3 = matrix(c(example_sd1^2,t_corr_init[1]*example_sd1*example_sd2,t_corr_init[2]*example_sd1*example_sd3,
                            t_corr_init[1]*example_sd1*example_sd2,example_sd2^2,t_corr_init[3]*example_sd2*example_sd3,
                            t_corr_init[2]*example_sd1*example_sd3,t_corr_init[3]*example_sd2*example_sd3,example_sd3^2),
                          3,3)
  
  
      data_example2 = data.frame(signal = c(t(replicate(i_n,c(0,1,10)))) ,
                            rbind(
                              mvtnorm::rmvnorm(i_n,mean=example_means1,sigma=example_sigma1),
                              mvtnorm::rmvnorm(i_n,mean=example_means2,sigma=example_sigma2),
                              mvtnorm::rmvnorm(i_n,mean=example_means3,sigma=example_sigma3)
                            ))
      
      devtools::use_data(data_example2,overwrite=TRUE)
