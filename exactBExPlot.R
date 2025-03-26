source("https://raw.githubusercontent.com/icascos/expdepth/master/exactexp.R")

conf.level <- function(beta=.95,n=100) {
  sq.dist <- function(alpha,ee=-sqrt(qchisq(beta,df=2)/n)){(dnorm(ee)/ee+pnorm(ee)-alpha/(2*alpha-1))^2}
  return(optimize(f=sq.dist,interval=c(0,1))$minimum)
}


exactBExPlot <- function(data,conf=0.95) {
  plot(data,type="n",main="BExPlot",xlab="X",ylab="Y")
  av <- colMeans(data)
  bag <- exactexp(data,alpha=.15)
  alpha.conf <- conf.level(beta=conf,n=nrow(data))
  conf.region <- exactexp(data,alpha=alpha.conf)
  expanded <- t(4*(t(bag)-av)+av)
  polygon(bag,col="light grey",border=FALSE)
  lines(bag,lty=2)
  polygon(conf.region,col="dark grey",border=FALSE)
  lines(conf.region,lty=2)
  points(x=av[1],y=av[2],pch=19)
  out <- NULL
  for(i in 1:nrow(data)){
    if(max(chull(rbind(expanded,data[i,])))==(nrow(expanded)+1)){out<-c(out,i)}
  }
  if(!is.null(out)) {points(x=data[out,1],y=data[out,2]);data.in <- data[-out,]}
  if(is.null(out)){data.in <- data}
  points(data.in,pch=3,cex=.5)
  fence <- chull(data.in)
  fence <- c(fence,fence[1])
  lines(data.in[fence,])
}


# Some examples

#require(ddalpha)
#data(hemophilia)
#attach(hemophilia)
#AHF.normal <-cbind(AHFactivity[gr=="normal"],AHFactivity.1[gr=="normal"])
#exactBExPlot(AHF.normal)

#AHF.carrier <-cbind(AHFactivity[gr=="carrier"],AHFactivity.1[gr=="carrier"])
#exactBExPlot(AHF.carrier,conf=0.99)

# Run only if interested in testing the complexity of the algorithm

#data10 <- cbind(rnorm(10),rnorm(10))
#data50 <- cbind(rnorm(50),rnorm(50))
#data100 <- cbind(rnorm(100),rnorm(100))
#data500 <- cbind(rnorm(500),rnorm(500))
#data1000 <- cbind(rnorm(1000),rnorm(1000))
#require(microbenchmark)
#microbenchmark(exactexp(data=data10),
#               exactexp(data=data50),
#               exactexp(data=data100),
#               exactexp(data=data500),times=20)
