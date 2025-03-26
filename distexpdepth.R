# Some distortion functions
identity <- function(x){return(x)}
square <- function(x){return(x^2)}
sroot <- function(x){return(sqrt(x))}
trimf = function(x,beta=.1) {return(punif(x,min=beta/2,max=1-beta/2))}
sigmoid <- function(x,delta=3){return((x!=0)*(1/(1+(1/x-1)^(delta))))}


# Distorted expectile depth of x wrt data with dual distortion function gtilde 
distexpdepth <- function(x,data,gtilde, sim=F){
  # set sim=TURE if distortion function is symmetric, g=gtilde
  if(min(chull(rbind(x,data)))==1) {return(0)}
  
  n <- nrow(data)
  data <- t(t(data)-x)
  
  steps <- seq(from=1,to=0,length.out=n+1)
  w <- gtilde(steps[-(n+1)])-gtilde(steps[-1])
  
  gamma <- atan2(data[,2],data[,1])
  ANG1 <- cbind(1:n,1:n,gamma+pi*((gamma<(-pi/2))-(gamma>pi/2)))
  side <- 2*(gamma>=-pi/2)*(gamma<=pi/2)-1
  ANG2 <- matrix(nrow=choose(n,2),ncol=3)
  k <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      gamma <- atan2(data[j,2]-data[i,2],data[j,1]-data[i,1])
      ANG2[k,] <- c(i,j,gamma+pi*((gamma<(-pi/2))-(gamma>pi/2)))
      k <- k+1
    }
  }
  ANG <- rbind(ANG1,ANG2)
  ord <- order(ANG[,3])
  ANG <- ANG[ord,]
  
  # if g is different from gtile, we must consider the whole of te circumference
  if(sim==FALSE) ANG <- rbind(ANG,cbind(ANG[,1:2],ANG[,3]+pi))
  
  depth <- 0
  rel.pos <- rank(data[,1])
  
  av.sum.pos <- colSums(w[rel.pos]*data[,1:2]*(side==1))
  
  av.sum.neg <- colSums(w[rel.pos]*data[,1:2]*(side==-1))
  
  for(i in 1:nrow(ANG)) {
    depth1 <- -(av.sum.neg[1]*cos(ANG[i,3]+pi/2)+av.sum.neg[2]*sin(ANG[i,3]+pi/2))/(av.sum.pos[1]*cos(ANG[i,3]+pi/2)+av.sum.pos[2]*sin(ANG[i,3]+pi/2))
    if(ANG[i,1]==ANG[i,2]) {
      k <- ANG[i,1]
      av.sum.pos <- av.sum.pos+w[rel.pos[k]]*data[k,1:2]*(side[k]==-1)
      av.sum.neg <- av.sum.neg-w[rel.pos[k]]*data[k,1:2]*(side[k]==-1)
      av.sum.pos <- av.sum.pos-w[rel.pos[k]]*data[k,1:2]*(side[k]==1)
      av.sum.neg <- av.sum.neg+w[rel.pos[k]]*data[k,1:2]*(side[k]==1)
      side[k] <- -side[k]
    }
    if(ANG[i,1]!=ANG[i,2]) {
      if(side[ANG[i,1]]==1){
        av.sum.pos <- av.sum.pos-w[rel.pos[ANG[i,1]]]*data[ANG[i,1],1:2] 
        av.sum.pos <- av.sum.pos-w[rel.pos[ANG[i,2]]]*data[ANG[i,2],1:2] 
        aux <- rel.pos[ANG[i,1]]
        rel.pos[ANG[i,1]] <- rel.pos[ANG[i,2]]
        rel.pos[ANG[i,2]] <- aux
        av.sum.pos <- av.sum.pos+w[rel.pos[ANG[i,1]]]*data[ANG[i,1],1:2] 
        av.sum.pos <- av.sum.pos+w[rel.pos[ANG[i,2]]]*data[ANG[i,2],1:2]
      }
      if(side[ANG[i,1]]==-1){
        av.sum.neg <- av.sum.neg-w[rel.pos[ANG[i,1]]]*data[ANG[i,1],1:2] 
        av.sum.neg <- av.sum.neg-w[rel.pos[ANG[i,2]]]*data[ANG[i,2],1:2] 
        aux <- rel.pos[ANG[i,1]]
        rel.pos[ANG[i,1]] <- rel.pos[ANG[i,2]]
        rel.pos[ANG[i,2]] <- aux
        av.sum.neg <- av.sum.neg+w[rel.pos[ANG[i,1]]]*data[ANG[i,1],1:2] 
        av.sum.neg <- av.sum.neg+w[rel.pos[ANG[i,2]]]*data[ANG[i,2],1:2]
      }
    }

    depth2 <- -(av.sum.neg[1]*cos(ANG[i,3]+pi/2)+av.sum.neg[2]*sin(ANG[i,3]+pi/2))/(av.sum.pos[1]*cos(ANG[i,3]+pi/2)+av.sum.pos[2]*sin(ANG[i,3]+pi/2))
    depth <- max(depth,depth1,depth2)
    
    # if g=gtilde, a semicircumference is enough
    if(sim==TRUE) 
      depth <- max(depth,depth1,depth2,1/depth1,1/depth2)
  }  
  
  distexpdepth <- 1/(1+depth)
  return(distexpdepth)
}
