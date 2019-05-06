aux_loadPackages<-function(packs){
  
  for (ipack in packs){
    if(require(ipack,character.only = TRUE)){
      print(paste0(ipack,": loaded correctly"))
    } else {
      print(paste0("trying to install ",ipack))
      install.packages(ipack)
      if(require(ipack,character.only = TRUE)){
        print(paste0(ipack, ": installed and loaded"))
      } else {
        stop(paste0("could not install ", ipack))
      }
    }
  }
  
}

cc_exact<-function(type="gauss",classes=2,dimension=1,...) {
  
  if (type=="gauss") {
    if (classes==2){
      if (dimension==1){
        cc_exact_gauss_2Cl1D(...)
      } else if (dimension==2) {
        cc_exact_gauss_2Cl2D(...)
      } else {
        cc_exact_gauss_2ClmD(...)
      }
    } else if (classes==3){
      if (dimension==1){
        cc_exact_gauss_3Cl1D(...)
      } else if (dimension==3) {
        cc_exact_gauss_3Cl3D(...)
      } else {
        stop("wrong dimension given")
      }
    } else {
      if (dimension==1){
        cc_exact_gauss_kCl1D(...)
      } else {
        cc_exact_gauss_kClmD(...)
      }
    }
  } else if (type=="exp") {
    if (classes==2) {
    cc_exact_exp_2Cl(...)
    } else if (classes==3){
    cc_exact_exp_3Cl(...)
    } else {
    cc_exact_exp_kCl(...)
    }
  } else if (type=="gamma") {
    cc_exact_gamma_kCl(...)
  } else if (type=="lnorm") {
    cc_exact_lnorm_kCl(...)
  } else {
    stop('Wrong type of distribution')
  }
  
  
}



cc_exact_gauss_2Cl1D<-function(par,x_length,s_length){
  mu1=par$mu1
  mu2=par$mu2
  sd1=par$sd1
  sd2=par$sd2
  
  x_span=seq(from=(min(mu1,mu2)-5*max(sd1,sd2) ), to=(max(mu1,mu2)+5*max(sd1,sd2) ), length=x_length )
  s_span=seq(from=0,to=1,length=s_length)
  
  MI=c()
  
  for (s in 1:length(s_span)){
    
    px    = s_span[s]*dnorm(x_span,mean=mu1,sd=sd1) + (1-s_span[s])*dnorm(x_span,mean=mu2,sd=sd2)
    A     =    s_span[s]  * trapz_simple(x_span, log( px^dnorm(x_span,mean=mu1,sd=sd1) ) )
    B     = (1-s_span[s]) * trapz_simple(x_span, log( px^dnorm(x_span,mean=mu2,sd=sd2) ) )
    MI[s] = - ( s_span[s]*0.5*log(2*pi*exp(1)*sd1^2) + (1-s_span[s])*0.5*log(2*pi*exp(1)*sd2^2) ) - (A+B)
    
  }
  
  #plot(MI)
  #which.max(MI)
  output=list()
  output$s_span=s_span
  output$MI=log2(exp(MI))
  output$cc=log2(max(exp(MI),na.rm=TRUE))
  output
}



cc_exact_gauss_2Cl2D<-function(par,x_length,s_length){
  mu1=par$mu1
  mu2=par$mu2
  sigma1=par$sigma1
  sigma2=par$sigma2
  
  
  x1_span=seq(from=(min(mu1[1],mu2[1])-5*sqrt(max(sigma1[1,1],sigma2[1,1])) ),
              to=(max(mu1[1],mu2[1])+5*sqrt(max(sigma1[1,1],sigma2[1,1])) ), length=x_length )
  x2_span=seq(from=(min(mu1[2],mu2[2])-5*sqrt(max(sigma1[2,2],sigma2[2,2])) ),
              to=(max(mu1[2],mu2[2])+5*sqrt(max(sigma1[2,2],sigma2[2,2])) ), length=x_length )
  s_span=seq(from=0,to=1,length=s_length)
  
  MI=c()
  
  for (s in 1:length(s_span)){
    px=matrix(0,length(x1_span),length(x2_span))
    pxsx1=matrix(0,length(x1_span),length(x2_span))
    pxsx2=matrix(0,length(x1_span),length(x2_span))
    
    for (i in 1:length(x1_span)){
      for (j in 1:length(x2_span)){
        px[i,j]    = s_span[s]*(mvtnorm::dmvnorm(x=c(x1_span[i],x2_span[j]),mean=mu1,sigma=sigma1)) + 
          (1-s_span[s])*(mvtnorm::dmvnorm(x=c(x1_span[i],x2_span[j]),mean=mu2,sigma=sigma2))
        pxsx1[i,j]  = log(px[i,j] ^ (mvtnorm::dmvnorm(x=c(x1_span[i],x2_span[j]),mean=mu1,sigma=sigma1)) )
        pxsx2[i,j]  = log(px[i,j] ^ (mvtnorm::dmvnorm(x=c(x1_span[i],x2_span[j]),mean=mu2,sigma=sigma2)) )
      }
    }
    
    A     =    s_span[s]  * trapz_simple(x1_span, apply(pxsx1,1,function(xx) trapz_simple(x2_span,xx)  ))
    B     = (1-s_span[s]) * trapz_simple(x1_span, apply(pxsx2,1,function(xx) trapz_simple(x2_span,xx)  ))
    E     = - ( s_span[s]*0.5*log( (2*pi*exp(1))^2* det(sigma1) ) + (1-s_span[s])*0.5*log( (2*pi*exp(1))^2 *det(sigma2)  ) )
    MI[s] = E - (A+B)
    
  }
  
  #plot(MI)
  #which.max(MI)
  output=list()
  output$s_span=s_span
  output$MI=log2(exp(MI))
  output$cc=log2(max(exp(MI),na.rm=TRUE))
  output
}



cc_exact_gauss_3Cl1D<-function(par,x_length,s_length){
  mu1=par$mu1
  mu2=par$mu2
  mu3=par$mu3
  sd1=par$sd1
  sd2=par$sd2
  sd3=par$sd3
  
  s_span=seq(from=0,to=1,length=s_length)
  MI=matrix(NA,length(s_span),length(s_span))  
  
  x_span=seq(from=(min(mu1,mu2,mu3)-4*max(sd1,sd2,sd3) ), to=(max(mu1,mu2,mu3)+4*max(sd1,sd2,sd3) ), length=x_length )
  
  for (s1 in 1:length(s_span)){
    for (s2 in 1:length(s_span)){    
      if (s_span[s1]+s_span[s2]< 1){
        px    = s_span[s1]*dnorm(x_span,mean=mu1,sd=sd1) + s_span[s2]*dnorm(x_span,mean=mu2,sd=sd2)+ (1-s_span[s1]-s_span[s2])*dnorm(x_span,mean=mu3,sd=sd3)
        A     =    s_span[s1]  * trapz_simple(x_span, log( px^dnorm(x_span,mean=mu1,sd=sd1) ) )
        B     =    s_span[s2]  * trapz_simple(x_span, log( px^dnorm(x_span,mean=mu2,sd=sd2) ) )
        C     = (1-s_span[s1]-s_span[s2]) * trapz_simple(x_span, log( px^dnorm(x_span,mean=mu3,sd=sd3) ) )
        MI[s1,s2] = - ( s_span[s1]*0.5*log(2*pi*exp(1)*sd1^2) + s_span[s2]*0.5*log(2*pi*exp(1)*sd2^2) + (1-s_span[s1]-s_span[s2])*0.5*log(2*pi*exp(1)*sd3^2) ) - (A+B+C)
      } else {
        MI[s1,s2]= 0
      }
    }
  }
  
  output=list()
  output$s_span=s_span
  output$MI=log2(exp(MI))
  output$cc=log2(max(exp(MI),na.rm=TRUE))
  output
  #   which_max_j=which.max(apply(MI,2,max))
  #   which_max_i=apply(MI,2,which.max)[which_max_j]
}



cc_exact_gauss_3Cl3D<-function(par,x_length,s_length){
  mu1=par$mu1
  mu2=par$mu2
  mu3=par$mu3
  sigma1=par$sigma1
  sigma2=par$sigma2
  sigma3=par$sigma3
  
  
  x1_span=seq(from=(min(mu1[1],mu2[1],mu3[1])-5*sqrt(max(sigma1[1,1],sigma2[1,1],sigma3[1,1])) ),
              to=(max(mu1[1],mu2[1],mu3[1])+5*sqrt(max(sigma1[1,1],sigma2[1,1],sigma3[1,1])) ), length=x_length )
  x2_span=seq(from=(min(mu1[2],mu2[2],mu3[2])-5*sqrt(max(sigma1[2,2],sigma2[2,2],sigma3[2,2])) ),
              to=(max(mu1[2],mu2[2],mu3[2])+5*sqrt(max(sigma1[2,2],sigma2[2,2],sigma3[2,2])) ), length=x_length )
  x3_span=seq(from=(min(mu1[3],mu2[3],mu3[3])-5*sqrt(max(sigma1[3,3],sigma2[3,3],sigma3[3,3])) ),
              to=(max(mu1[3],mu2[3],mu3[3])+5*sqrt(max(sigma1[3,3],sigma2[3,3],sigma3[3,3])) ), length=x_length )
  s_span=seq(from=0,to=1,length=s_length)
  
  MI=matrix(NA,length(s_span),length(s_span))  
  
  for (s1 in 1:length(s_span)){
    print(paste("Now:",s1))
    for (s2 in 1:length(s_span)){
      if (s_span[s1]+s_span[s2]< 1){
        #         px=array(0,dim=c(length(x1_span),length(x2_span),length(x3_span)))
        #         pxsx1=array(0,dim=c(length(x1_span),length(x2_span),length(x3_span)))
        #         pxsx2=array(0,dim=c(length(x1_span),length(x2_span),length(x3_span)))
        #         pxsx3=array(0,dim=c(length(x1_span),length(x2_span),length(x3_span)))
        #         
        #         for (i in 1:length(x1_span)){
        #           for (j in 1:length(x2_span)){
        #             for (k in 1:length(x3_span)){
        #             px[i,j,k]    = s_span[s1]*dmvnorm(x=c(x1_span[i],x2_span[j],x3_span[k]),mean=mu1,sigma=sigma1) +
        #                            s_span[s2]*dmvnorm(x=c(x1_span[i],x2_span[j],x3_span[k]),mean=mu2,sigma=sigma2) + 
        #                            (1-s_span[s1]-s_span[s2])*dmvnorm(x=c(x1_span[i],x2_span[j],x3_span[k]),mean=mu3,sigma=sigma3)
        #             pxsx1[i,j,k]  = log(px[i,j,k] ^ dmvnorm(x=c(x1_span[i],x2_span[j],x3_span[k]),mean=mu1,sigma=sigma1) )
        #             pxsx2[i,j,k]  = log(px[i,j,k] ^ dmvnorm(x=c(x1_span[i],x2_span[j],x3_span[k]),mean=mu2,sigma=sigma2) )
        #             pxsx3[i,j,k]  = log(px[i,j,k] ^ dmvnorm(x=c(x1_span[i],x2_span[j],x3_span[k]),mean=mu3,sigma=sigma3) )
        #             }
        #           }
        #         }
        
        temp_X  = expand.grid(x1_span,x2_span,x3_span)
        temp_PX = s_span[s1]*(mvtnorm::dmvnorm(x=temp_X,mean=mu1,sigma=sigma1)) +
          s_span[s2]*(mvtnorm::dmvnorm(x=temp_X,mean=mu2,sigma=sigma2)) + 
          (1-s_span[s1]-s_span[s2])*(mvtnorm::dmvnorm(x=temp_X,mean=mu3,sigma=sigma3))
        pxsx1  = array(log(temp_PX ^ (mvtnorm::dmvnorm(x=temp_X,mean=mu1,sigma=sigma1)) ),
                       dim=c(x_length,x_length,x_length))
        pxsx2  = array(log(temp_PX ^ (mvtnorm::dmvnorm(x=temp_X,mean=mu2,sigma=sigma2)) ),
                       dim=c(x_length,x_length,x_length))
        pxsx3  = array(log(temp_PX ^ (mvtnorm::dmvnorm(x=temp_X,mean=mu3,sigma=sigma3)) ),
                       dim=c(x_length,x_length,x_length))
        
        
        A     =    s_span[s1]  * trapz_simple(x1_span, apply(  apply(pxsx1,1:2,function(xxx) trapz_simple(x3_span,xxx))  ,1,function(xx) trapz_simple(x2_span,xx)  ))
        B     =   (s_span[s2]) * trapz_simple(x1_span, apply(  apply(pxsx2,1:2,function(xxx) trapz_simple(x3_span,xxx))  ,1,function(xx) trapz_simple(x2_span,xx)  ))
        C     = (1-s_span[s1]-s_span[s2]) * trapz_simple(x1_span, apply(  apply(pxsx3,1:2,function(xxx) trapz_simple(x3_span,xxx))  ,1,function(xx) trapz_simple(x2_span,xx)  ))
        E     = - ( s_span[s1]*0.5*log( (2*pi*exp(1))^3* det(sigma1) ) + 
                      s_span[s2]*0.5*log( (2*pi*exp(1))^3* det(sigma2) ) + 
                      (1-s_span[s1]-s_span[s2])*0.5*log( (2*pi*exp(1))^3 *det(sigma3)  ) )
        MI[s1,s2] = E - (A+B+C)
      } else {
        MI[s1,s2]= 0
      }
    }
  }
  
  output=list()
  output$s_span=s_span
  output$MI=log2(exp(MI))
  output$cc=log2(max(exp(MI),na.rm=TRUE))
  output
  #   which_max_j=which.max(apply(MI,2,max))
  #   which_max_i=apply(MI,2,which.max)[which_max_j]
}



cc_exact_gauss_kCl1D<-function(par,x_length,mc_length,mc_cores=1){
  mu=par$mu
  sd=par$sd
  
  `%dopar%`<-foreach::`%dopar%`
  MI=c()
  s_samp=c()
  k=length(mu)
  x_span=seq(from=(min(mu)-10*max(sd) ), to=(max(mu)+10*max(sd) ), length=x_length )
  # N_samp=rep(N,k)
  
  # dataTest  <- data.frame( response=do.call(c,mapply( function(x,y,z){rnorm(z,mean=x,sd=y)},mu,sd,N_samp,SIMPLIFY=FALSE)), signal=do.call(c,mapply(function(x,y) rep(x,y),letters[1:k],N_samp,SIMPLIFY=FALSE )) )
  #   xklr=dataTest[,1]
  #   yklr=dataTest[,2]
  #   ggplot(data=dataTest,aes(x=response,y=..density..))+geom_histogram(aes(colour=signal))+facet_grid(signal~.)
  #   
  
  #dataMI<-data.frame(response=xklr,signal=yklr)
  
  cl=parallel::makeCluster(mc_cores)
  doParallel::registerDoParallel(cl)
  OutList=foreach::foreach(i=1:mc_length,.export=c("trapz_simple","x_log_y")) %dopar% {
    out_temp<-list()
    s      =   runif(k)
    s      =   s/sum(s)
    px     =    apply(mapply(function(x,y,z) {z*dnorm(x_span,mean=x,sd=y)} ,mu,sd,s,SIMPLIFY=TRUE ),1,sum)
    Ilogpx =   mapply(function(x,y,z) {z* trapz_simple(x_span, x_log_y( dnorm(x_span,mean=x,sd=y),px ) ) } ,mu,sd,s,SIMPLIFY=TRUE )
    
    out_temp$MI = - sum(s*0.5*log(2*pi*exp(1)*sd^2))  - sum(Ilogpx)
    out_temp$s_samp=s
    out_temp
  }
  parallel::stopCluster(cl)
  
  output=list()
  output$s_span=t(sapply(OutList,function(x) x$s_samp))
  output$MI=sapply(OutList,function(x) x$MI)
  output$cc=log2(max(exp(output$MI),na.rm=TRUE))
  output
}



cc_exact_gauss_2ClmD<-function(par,s_length,mc2_length,mc_cores=1){
  mu1=par$mu1
  sigma1=par$sigma1
  mu2=par$mu2
  sigma2=par$sigma2
  
  mu=list()
  mu[[1]]=mu1
  mu[[2]]=mu2
  
  sigma=list()
  sigma[[1]]=sigma1
  sigma[[2]]=sigma2
  
  `%dopar%`<-foreach::`%dopar%`
  MI=c()
  s_span=seq(from=0,to=1,length=s_length)
  m=length(mu[[1]])
  
  cl=parallel::makeCluster(mc_cores)
  doParallel::registerDoParallel(cl)
  OutList=foreach::foreach(s=s_span,.packages = c("mvtnorm"),.export=c("trapz_simple","jointProb_gaussian_md","x_log_y")) %dopar% {
    temp_mcvalue=c()
    for (j in 1:mc2_length){
      temp_samp=sample(1:2,prob=c(s,1-s),size=1,replace=TRUE)
      temp_x=mvtnorm::rmvnorm(n=1,mean=mu[[temp_samp]],sigma=sigma[[temp_samp]])
      temp_mcvalue[j]=log(jointProb_gaussian_md(temp_x,c(s,1-s),mu,sigma))
    }
    Ilogpx =   mean(temp_mcvalue)
    out_temp<-list()
    out_temp$MI = - s*0.5*log( (2*pi*exp(1))^m * det(sigma1))- (1-s)*0.5*log( (2*pi*exp(1))^m * det(sigma2)) - Ilogpx
    out_temp$s_samp=s
    out_temp
  }
  parallel::stopCluster(cl)
  
  output=list()
  output$s_span=(sapply(OutList,function(x) x$s_samp))
  output$MI=sapply(OutList,function(x) x$MI)
  output$cc=log2(max(exp(output$MI),na.rm=TRUE))
  output
}




## Multi-Class
cc_exact_gauss_kClmD<-function(par,x_length,mc_length,mc2_length,mc_cores=1,popt_sugg=NULL){
  mu=par$mu
  sigma=par$sigma
  
  `%dopar%`<-foreach::`%dopar%`
  MI=c()
  s_samp=c()
  k=length(mu)
  m=length(mu[[1]])
  x_span=lapply(1:k,function(x) {
    seq(from=(min( sapply(mu,function(xx) xx[x])  )-5*max(sapply(sigma,function(xx) xx[x,x]) ) ), 
        to  =(max( sapply(mu,function(xx) xx[x]) ) +5*max(sapply(sigma,function(xx) xx[x,x]) ) ), 
        length=x_length )
  })
  # N_samp=rep(N,k)
  
  # dataTest  <- data.frame( response=do.call(c,mapply( function(x,y,z){rnorm(z,mean=x,sd=y)},mu,sd,N_samp,SIMPLIFY=FALSE)), signal=do.call(c,mapply(function(x,y) rep(x,y),letters[1:k],N_samp,SIMPLIFY=FALSE )) )
  #   xklr=dataTest[,1]
  #   yklr=dataTest[,2]
  #   ggplot(data=dataTest,aes(x=response,y=..density..))+geom_histogram(aes(colour=signal))+facet_grid(signal~.)
  #   
  
  #dataMI<-data.frame(response=xklr,signal=yklr)
  
  cl=parallel::makeCluster(mc_cores)
  doParallel::registerDoParallel(cl)
  OutList=foreach::foreach(i=1:mc_length,.packages = c("mvtnorm"),.export=c("trapz_simple","jointProb_gaussian_md","x_log_y")) %dopar% {
    out_temp<-list()
    s      =   runif(k)
    s      =   s/sum(s)
    temp_mcvalue=c()
    for (j in 1:mc2_length){
      temp_samp=sample(1:k,prob=s,size=1,replace=TRUE)
      temp_x=mvtnorm::rmvnorm(n=1,mean=mu[[temp_samp]],sigma=sigma[[temp_samp]])
      temp_mcvalue[j]=log(jointProb_gaussian_md(temp_x,s,mu,sigma))
    }
    Ilogpx =   mean(temp_mcvalue)
    out_temp$MI = - sum(s*0.5*log( (2*pi*exp(1))^m * sapply(sigma,function(x) det(x))  ))  - Ilogpx
    out_temp$s_samp=s
    out_temp
  }
  parallel::stopCluster(cl)
  
  output=list()
  output$s_span=t(sapply(OutList,function(x) x$s_samp))
  output$MI=sapply(OutList,function(x) x$MI)
  
  
  if (!is.null(popt_sugg)){
    s=popt_sugg/sum(popt_sugg)
    temp_mcvalue=c()
    for (j in 1:mc2_length){
      temp_samp=sample(1:k,prob=s,size=1,replace=TRUE)
      temp_x=mvtnorm::rmvnorm(n=1,mean=mu[[temp_samp]],sigma=sigma[[temp_samp]])
      temp_mcvalue[j]=log(jointProb_gaussian_md(temp_x,s,mu,sigma))
    }
    Ilogpx =   mean(temp_mcvalue)
    output$MI = c(output$MI,(- sum(s*0.5*log( (2*pi*exp(1))^m * sapply(sigma,function(x) det(x))  ))  - Ilogpx))
    output$s_span=rbind(output$s_span,s)
  }
  
  output$cc=log2(max(exp(output$MI),na.rm=TRUE))
  output
}



## EXPONTENTIALS

cc_exact_exp_2Cl<-function(par,x_length,s_length){
  lam1=par$lam1
  lam2=par$lam2
  
  
  MI=c()
  
  x_span=seq(from=0, to=3*log(10)*(1/min(lam1,lam2)), length=x_length)
  s_span=seq(from=0,to=1,length=s_length)
  
  for (s in 1:length(s_span)){
    
    px    = s_span[s]*dexp(x_span,rate=lam1) + (1-s_span[s])*dexp(x_span,rate=lam2)
    A     =    s_span[s]  * trapz_simple(x_span, log( px^dexp(x_span,rate=lam1) ) )
    B     = (1-s_span[s]) * trapz_simple(x_span, log( px^dexp(x_span,rate=lam2) ) )
    MI[s] = s_span[s]*log(lam1/lam2)+log(lam2)-1 - (A+B)
    
  }
  
  output=list()
  output$s_span=s_span
  output$MI=log2(exp(MI))
  output$cc=log2(max(exp(MI),na.rm=TRUE))
  output
}



cc_exact_exp_3Cl<-function(par,x_length,s_length){
  lam1=par$lam1
  lam2=par$lam2
  lam3=par$lam3
  
  s_span=seq(from=0,to=1,length=s_length)
  MI=matrix(NA,length(s_span),length(s_span))
  
  x_span=seq(from=0, to=3*log(10)*(1/min(lam1,lam2,lam3)), length=x_length)
  
  for (s1 in 1:length(s_span)){
    for (s2 in 1:length(s_span)){
      if (s_span[s1]+s_span[s2]< 1){
        px    = s_span[s1]*dexp(x_span,rate=lam1) + s_span[s2]*dexp(x_span,rate=lam2)+ (1- s_span[s1]- s_span[s2])*dexp(x_span,rate=lam3)
        A     =    s_span[s1]  * trapz_simple(x_span, log( px^dexp(x_span,rate=lam1) ) )
        B     =    s_span[s2]  * trapz_simple(x_span, log( px^dexp(x_span,rate=lam2) ) )
        C     = (1-s_span[s1]-s_span[s2]) * trapz_simple(x_span, log( px^dexp(x_span,rate=lam3) ) )
        MI[s1,s2] = s_span[s1]*(log(lam1)-1) + s_span[s2]*(log(lam2)-1) + (1-s_span[s1]-s_span[s2])*(log(lam3)-1)  - (A+B+C)
      } else {
        MI[s1,s2]= 0
      }
    }
  }
  
  filled.contour(s_span,s_span,MI)
  log2(max(exp(MI),na.rm=TRUE))
  # which_max_j=which.max(apply(MI,2,max))
  # which_max_i=apply(MI,2,which.max)[which_max_j]
}



cc_exact_exp_kCl<-function(par,x_length,mc_length,mc_cores=1){
  lam=par$lam
  
  MI=c()
  s_samp=c()
  k=length(lam)
  x_span=c(0,exp(seq(from=log(1e-6), to=log( qexp(1-1e-6,rate=min(lam)) ), length=x_length )))
  # N_samp=rep(N,k)
  
  # dataTest  <- data.frame( response=do.call(c,mapply( function(x,y,z){rnorm(z,mean=x,sd=y)},mu,sd,N_samp,SIMPLIFY=FALSE)), signal=do.call(c,mapply(function(x,y) rep(x,y),letters[1:k],N_samp,SIMPLIFY=FALSE )) )
  #   xklr=dataTest[,1]
  #   yklr=dataTest[,2]
  #   ggplot(data=dataTest,aes(x=response,y=..density..))+geom_histogram(aes(colour=signal))+facet_grid(signal~.)
  #   
  
  #dataMI<-data.frame(response=xklr,signal=yklr)
  
  `%dopar%`<-foreach::`%dopar%`
  
  cl=parallel::makeCluster(mc_cores)
  doParallel::registerDoParallel(cl)
  OutList=foreach::foreach(i=1:mc_length,.export=c("trapz_simple","x_log_y")) %dopar% {
    out_temp<-list()
    s      =   runif(k)
    s      =   s/sum(s)
    px     =    apply(mapply(function(x,z) {z*dexp(x_span,rate=x)} ,lam,s,SIMPLIFY=TRUE ),1,sum)
    Ilogpx =   mapply(function(x,z) {z* trapz_simple(x_span, log( px^dexp(x_span,rate=x) ) ) } ,lam,s,SIMPLIFY=TRUE )
    Ilogpxs =   mapply(function(x,z) {z* trapz_simple(x_span, log( dexp(x_span,rate=x)^dexp(x_span,rate=x) ) ) } ,lam,s,SIMPLIFY=TRUE )
    
    out_temp$MI = sum(Ilogpxs)  - sum(Ilogpx)
    out_temp$s_samp=s
    out_temp
  }
  parallel::stopCluster(cl)
  
  output=list()
  output$s_span=t(sapply(OutList,function(x) x$s_samp))
  output$MI=sapply(OutList,function(x) x$MI)
  output$cc=log2(max(exp(output$MI),na.rm=TRUE))
  output
}





cc_exact_lnorm_kCl<-function(par,x_length,mc_length,mc_cores=1){
  mu=par$mu
  sd=par$sd
  
  MI=c()
  s_samp=c()
  k=length(mu)
  x_span=exp(seq(from=min(mu-5*sd),to=max(mu+5*sd),length.out = x_length))
  # N_samp=rep(N,k)
  
  # dataTest  <- data.frame( response=do.call(c,mapply( function(x,y,z){rnorm(z,mean=x,sd=y)},mu,sd,N_samp,SIMPLIFY=FALSE)), signal=do.call(c,mapply(function(x,y) rep(x,y),letters[1:k],N_samp,SIMPLIFY=FALSE )) )
  #   xklr=dataTest[,1]
  #   yklr=dataTest[,2]
  #   ggplot(data=dataTest,aes(x=response,y=..density..))+geom_histogram(aes(colour=signal))+facet_grid(signal~.)
  #   
  
  #dataMI<-data.frame(response=xklr,signal=yklr)
  
  `%dopar%`<-foreach::`%dopar%`
  
  cl=parallel::makeCluster(mc_cores)
  doParallel::registerDoParallel(cl)
  OutList=foreach::foreach(i=1:mc_length,.export=c("trapz_simple","x_log_y")) %dopar% {
    out_temp<-list()
    s      =   runif(k)
    s      =   s/sum(s)
    px     =    apply(mapply(function(x,y,z) {z*dlnorm(x_span,meanlog =x,sdlog=y)} ,mu,sd,s,SIMPLIFY=TRUE ),1,sum)
    Ilogpx =   mapply(function(x,y,z) {z*trapz_simple(x_span, x_log_y(dlnorm(x_span,meanlog =x,sdlog=y),px) ) } ,mu,sd,s,SIMPLIFY=TRUE )
    Ilogpxs =   mapply(function(x,y,z) {z*trapz_simple(x_span, x_log_y(dlnorm(x_span,meanlog =x,sdlog=y),dlnorm(x_span,meanlog =x,sdlog=y) ) ) } ,mu,sd,s,SIMPLIFY=TRUE )
    
    out_temp$MI = sum(Ilogpxs)  - sum(Ilogpx)
    out_temp$s_samp=s
    out_temp
  }
  parallel::stopCluster(cl)
  
  output=list()
  output$s_span=t(sapply(OutList,function(x) x$s_samp))
  output$MI=sapply(OutList,function(x) x$MI)
  output$cc=log2(max(exp(output$MI),na.rm=TRUE))
  output
}


## GAMMA

cc_exact_gamma_kCl<-function(par,x_length,mc_length,mc_cores=1){
  shape=par$shape
  scale=par$scale
  
  
  MI=c()
  s_samp=c()
  k=length(shape)
  
  temp_minx=min(mapply(function(x,y) qgamma((1e-6),shape=x,scale=y),shape,scale,SIMPLIFY = TRUE))
  if(temp_minx==0){
    temp_minx=1e-150
  } else {
    temp_minx=temp_minx
  } 
  
  x_span=c(exp(seq(from=log(temp_minx),
                     to=log(max(mapply(function(x,y) qgamma((1-1e-6),shape=x,scale=y),shape,scale,SIMPLIFY = TRUE))),length.out=x_length)))
  # N_samp=rep(N,k)
  
  # dataTest  <- data.frame( response=do.call(c,mapply( function(x,y,z){rnorm(z,mean=x,sd=y)},mu,sd,N_samp,SIMPLIFY=FALSE)), signal=do.call(c,mapply(function(x,y) rep(x,y),letters[1:k],N_samp,SIMPLIFY=FALSE )) )
  #   xklr=dataTest[,1]
  #   yklr=dataTest[,2]
  #   ggplot(data=dataTest,aes(x=response,y=..density..))+geom_histogram(aes(colour=signal))+facet_grid(signal~.)
  #   
  
  #dataMI<-data.frame(response=xklr,signal=yklr)
  
  `%dopar%`<-foreach::`%dopar%`
  
  cl=parallel::makeCluster(mc_cores)
  doParallel::registerDoParallel(cl)
  OutList=foreach::foreach(i=1:mc_length,.export=c("trapz_simple","x_log_y")) %dopar% {
    out_temp<-list()
    s      =   runif(k)
    s      =   s/sum(s)
    px     =    apply(mapply(function(x,y,z) {z*dgamma(x_span,shape=x,scale=y)} ,shape,scale,s,SIMPLIFY=TRUE ),1,sum)
    Ilogpx =   mapply(function(x,y,z) {z* trapz_simple(x_span, x_log_y(dgamma(x_span,shape=x,scale=y),px ) ) } ,shape,scale,s,SIMPLIFY=TRUE )
    Ilogpxs =   mapply(function(x,y,z) {z* trapz_simple(x_span, x_log_y(dgamma(x_span,shape=x,scale=y),dgamma(x_span,shape=x,scale=y) ) ) } ,shape,scale,s,SIMPLIFY=TRUE )
    
    out_temp$MI = sum(Ilogpxs)  - sum(Ilogpx)
    out_temp$s_samp=s
    out_temp
  }
  parallel::stopCluster(cl)
  
  output=list()
  output$s_span=t(sapply(OutList,function(x) x$s_samp))
  output$MI=sapply(OutList,function(x) x$MI)
  output$cc=log2(max(exp(output$MI),na.rm=TRUE))
  output
}


x_log_y<-function(x,y){
  out=log(y^x)
  ids=is.infinite(out)
  out[ids]=x[ids]*log(y[ids])
  out
}

trapz_simple<-function(x,y){
  nx=length(x)
  
  xmesh=x[(2:nx)]-x[(1:(nx-1))]
  ymesh=0.5*(y[(2:nx)]+y[(1:(nx-1))])
  
  sum(xmesh*ymesh)
}

jointProb_gaussian_md<-function(x,s,mu,sigma){
  out<-do.call(sum,mapply(function(signal,tmp_mean,tmp_sigma) {signal*mvtnorm::dmvnorm(x=x,mean=tmp_mean,sigma=tmp_sigma)},
                          s,mu,sigma,SIMPLIFY=FALSE))
  
  out
}



tempFun2_plot<-function(OutputList,i_class,i_output,path_output_main,plot_height,plot_width,ylimits=c(0,0.5),sdlimit=5){
  outputDF<-melt(lapply(OutputList,function(x) {
    temp_capacityexact=x$exact$cc
    lapply(x,function(xx) {
      if (is.null(xx$cc)){
        tempcc=(sapply(xx,function(xxx) abs(xxx$cc-temp_capacityexact)/temp_capacityexact ))
        data.frame(capacity=mean(tempcc),
                   error=sd(tempcc)  )
      } else {
        data.frame(capacity=0,error=0 )
      }
    }) 
  }),id.vars=c("capacity","error"))
  colnames(outputDF)<-c("capacity","error","N","sd")
  outputDF$sd<-as.numeric(outputDF$sd)
  
  
  PlotCapacityComparison=ggplot(data=outputDF[outputDF$N%in%c("100","500","1000")&(outputDF$sd<=sdlimit),],aes(x=sd,y=capacity,colour=N))+
    geom_point(size=1)+geom_line(size=0.5,alpha=0.5)+
    geom_errorbar(aes(ymin = capacity-error, ymax = capacity+error), width = 0.1,size=0.5)+
    scale_y_continuous("Relative error (%)",limits=ylimits)+scale_x_continuous(paste("Stand Dev.",sep="") )+
    ggtitle(paste("Capacity error -","Class:",i_class, "Dim:",i_output,sep=" " ) )+
    theme_publ(version=2)
  
  
  filepath=paste(path_output_main,"plotComparison.pdf",sep="")
  ggsave(PlotCapacityComparison,file=filepath,height=plot_height,width=plot_width)
  filepath=paste(path_output_main,"plotComparison_logx.pdf",sep="")
  ggsave(PlotCapacityComparison+scale_x_log10(),file=filepath,height=plot_height,width=plot_width)
  filepath=paste(path_output_main,"plotComparison_logxy.pdf",sep="")
  ggsave(PlotCapacityComparison+scale_x_log10()+scale_y_log10(),file=filepath,height=plot_height,width=plot_width)
  
  PlotCapacityComparison+scale_x_log10()
}

tempFun3_plot<-function(OutputList,i_class,i_output,path_output_main,plot_height,plot_width,ylimits=c(0,0.5),sdlimit=5){
  outputDF<-melt(lapply(OutputList,function(x) {
    temp_capacityexact=x$exact$cc
    lapply(x,function(xx) {
      if (is.null(xx$cc)){
        tempcc=(sapply(xx,function(xxx) abs(xxx$cc-temp_capacityexact) ))
        data.frame(capacity=mean(tempcc),
                   error=sd(tempcc)  )
      } else {
        data.frame(capacity=0,error=0 )
      }
    }) 
  }),id.vars=c("capacity","error"))
  colnames(outputDF)<-c("capacity","error","N","sd")
  outputDF$sd<-as.numeric(outputDF$sd)
  
  
  PlotCapacityComparison=ggplot(data=outputDF[outputDF$N%in%c("100","500","1000")&(outputDF$sd<=sdlimit),],aes(x=sd,y=capacity,colour=N))+
    geom_point(size=1)+geom_line(size=0.5,alpha=0.5)+
    geom_errorbar(aes(ymin = capacity-error, ymax = capacity+error), width = 0.1,size=0.5)+
    scale_y_continuous("Absolute error (bits)",limits=ylimits)+scale_x_continuous(paste("Stand Dev.",sep="") )+
    ggtitle(paste("Capacity error -","Class:",i_class, "Dim:",i_output,sep=" " ) )+
    theme_publ(version=2)

  
    filepath=paste(path_output_main,"plotAbComparison.pdf",sep="")
  ggsave(PlotCapacityComparison,file=filepath,height=plot_height,width=plot_width)
  filepath=paste(path_output_main,"plotAbComparison_logx.pdf",sep="")
  ggsave(PlotCapacityComparison+scale_x_log10(),file=filepath,height=plot_height,width=plot_width)
  filepath=paste(path_output_main,"plotAbComparison_logxy.pdf",sep="")
  ggsave(PlotCapacityComparison+scale_x_log10()+scale_y_log10(),file=filepath,height=plot_height,width=plot_width)
  
  PlotCapacityComparison+scale_x_log10()
}

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
