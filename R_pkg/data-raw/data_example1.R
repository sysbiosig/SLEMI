signal=c("0","0.1","1")
means=list()
sds=list()
means[["0"]]=   c(1, 2  ,  3 )
means[["0.1"]]= c(2, 4  ,  6 )
means[["1"]]=   c(3, 6  , 9   )

sds[["0"]]   = c(1   , 1 ,1)
sds[["0.1"]] = c(1 , 1   ,1  )
sds[["1"]]   = c(1 ,1 ,1)

N=c(500,500,500)
names(N)<-signal

sideTypes=c("type 1","type 2","type 3")
probs=c(0.6,0.3,0.1)

set.seed(123)
data_example1=do.call(rbind,lapply(signal,function(x){
  components<-sample(1:3,prob=probs,size=N[x],replace=TRUE)
  sideVarSample <- sideTypes[components]
  ResponseSample <- rnorm(n=N[x],mean=means[[x]][components],sd=sds[[x]][components])
  data.frame(sideVar=sideVarSample,response=ResponseSample,signal=x)
}))

colnames(data_example1)<-c("sideVar","response","signal")

devtools::use_data(data_example1,overwrite=TRUE)