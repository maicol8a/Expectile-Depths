source("https://raw.githubusercontent.com/cran/expectreg/master/R/expectile.R")


extreme.points <- function(data,nlines=1000,alpha=0.15){
  angle <- seq(0,2*pi,length.out=nlines) 
  extreme.points <- NULL
  exp1 <- expectile(x=data[,1]*cos(angle[1])+data[,2]*sin(angle[1]),probs=alpha,dec=16)
  for(i in 2:nlines){
    exp2 <-  expectile(x=data[,1]*cos(angle[i])+data[,2]*sin(angle[i]),probs=alpha,dec=16)
    extreme.points <- rbind(extreme.points,solve(matrix(c(cos(angle[i-1]),sin(angle[i-1]),
                                                          cos(angle[i]),sin(angle[i])),byrow=TRUE,ncol=2),
                                                 c(exp1,exp2)))
    exp1 <- exp2
  }
  extreme.points <- extreme.points[chull(extreme.points),]
  extreme.points <- rbind(extreme.points,extreme.points[1,])
  return(extreme.points)
}


conf.level <- function(beta=.95,n=100) {
  sq.dist <- function(alpha,ee=-sqrt(qchisq(beta,df=2)/n)){(dnorm(ee)/ee+pnorm(ee)-alpha/(2*alpha-1))^2}
  return(optimize(f=sq.dist,interval=c(0,1))$minimum)
}


BExPlot <- function(data,nlines=1000,conf=0.95) {
  plot(data,type="n",main="BExPlot",xlab="X",ylab="Y")
  av <- colMeans(data)
  bag <- extreme.points(data,nlines=nlines)
  alpha.conf <- conf.level(beta=conf,n=nrow(data))
  conf.region <- extreme.points(data,nlines=nlines,alpha=alpha.conf)
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
