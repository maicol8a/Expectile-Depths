source("https://raw.githubusercontent.com/cran/expectreg/master/R/expectile.R")

################################################

exactexp <- function(data,alpha=.15){
  # Step 1
  y  <- data[order(data[,1]),2]
  x <- sort(data[,1])
  # Step 4
  expect <- expectile(x,probs=alpha,dec=16)
  
  U <- upper_boundary(x,y,alpha,expect)
  L <- upper_boundary(x,-y,alpha,expect)
  
  L <- L%*%matrix(c(1,0,0,-1),ncol=2)
  E <- rbind(U,L[nrow(L):1,])
  return(E)
}

################################################

upper_boundary <- function(x,y,alpha,expect){
  n <- length(x)
  mx <- mean(x)
  my <- mean(y)
  # Step 2
  R<-(1:n)
  iR<-(1:n)
  # Step 3
  ANG <- xy.angles(x,y)
  s<-sum(x<=expect)  
#  Step 5
  EXT<-NULL

i<-1
M<-NULL
#Determine the candidates to extreme points 
s1<-s
s2<-s

# Step 6
#Extreme point for s
#1
y0<-vector_mas(s1,x,y,mx,my,alpha,iR)
# Step 7
M<-rbind(M,y0)

#Candidate for s+1
s1<-s+1
#2

# Step 9
y1<-vector_masp(s1,x,y,mx,my,alpha,iR,y0)
while (check.point(s1,x,y,0,ANG[i,1],y1,iR)) {
  M<-rbind(M,y1)
  s1 <- s1+1
  if(s1>=n)break()
  #3
  y1 <- vector_masp(s1,x,y,mx,my,alpha,iR,y1)
}

# Step 9
#Candidate for s-1
t1 <- s-1
#4
y2<-vector_masm(t1,x,y,mx,my,alpha,iR,y0)
while(t1>0){
  if(!check.point(t1,x,y,0,ANG[i,1],y2,iR))break()
  M<-rbind(M,y2)
  t1 <- t1-1
  #5
  y2<-vector_masm(t1,x,y,mx,my,alpha,iR,y2)
}


# Step 10
if(!is.null(M)) M = matrix(M[order(M[,4], decreasing = F),], ncol = 6)
s=ifelse(is.null(M), s, M[nrow(M),3])


y0=M[nrow(M),]

if((R[ANG[i,2]]>=(s+1))&(R[ANG[i,3]]<(s+1))) {
  y0[5]=y0[5]-x[ANG[i,2]]+x[ANG[i,3]]
  y0[6]=y0[6]-y[ANG[i,2]]+y[ANG[i,3]]
    }
if((R[ANG[i,3]]>=(s+1))&(R[ANG[i,2]]<(s+1))) {
  y0[5]=y0[5]-x[ANG[i,3]]+x[ANG[i,2]]
  y0[6]=y0[6]-y[ANG[i,3]]+y[ANG[i,2]]
}

# Step 11
aux=R[ANG[i,2]]
R[ANG[i,2]]=R[ANG[i,3]]
R[ANG[i,3]]=aux

iR[R[ANG[i,2]]]=ANG[i,2]
iR[R[ANG[i,3]]]=ANG[i,3]



EXT=rbind(EXT,M[,1:2])

for(i in 2:((n*(n-1)/2))) {
  M=NULL
  # Step 7
  #6
  y0=vector_mase(s,x,y,mx,my,alpha,iR,y0)
  if((EXT[nrow(EXT),1]!=y0[1])|(EXT[nrow(EXT),2]!=y0[2])){
    # Step 8
    if(check.point(s,x,y,ANG[i-1,1],ANG[i,1],y0,iR)){
      M=rbind(M,y0)
    }
  }
  
  s1=s+1
  #7
  # Step 9
  y1=vector_masp(s1,x,y,mx,my,alpha,iR,y0)
  while (check.point(s1,x,y,ANG[i-1,1],ANG[i,1],y1,iR)) {
    M=rbind(M,y1)
    s1 = s1+1
    if(s1>=n)break()
    #8
    y1 = vector_masp(s1,x,y,mx,my,alpha,iR,y1)
  }
  
  t1 = s-1
  #9
  # Step 9
  y2=vector_masm(t1,x,y,mx,my,alpha,iR,y0)
  while(t1>0){
    if(!check.point(t1,x,y,ANG[i-1,1],ANG[i,1],y2,iR))break()
    M=rbind(M,y2)
    t1 = t1-1
    #10
    y2=vector_masm(t1,x,y,mx,my,alpha,iR,y2)
  }
  if(!is.null(M)) M = matrix(M[order(M[,4], decreasing = F),], ncol = 6)
  s=ifelse(is.null(M), s, M[nrow(M),3])

  # Step 10
  
  if(!is.null(M)){
    y0=M[nrow(M),]
  
    if((R[ANG[i,2]]>=(s+1))&(R[ANG[i,3]]<(s+1))) {
      y0[5]=y0[5]-x[ANG[i,2]]+x[ANG[i,3]]
      y0[6]=y0[6]-y[ANG[i,2]]+y[ANG[i,3]]
    }
    if((R[ANG[i,3]]>=(s+1))&(R[ANG[i,2]]<(s+1))) {
      y0[5]=y0[5]-x[ANG[i,3]]+x[ANG[i,2]]
      y0[6]=y0[6]-y[ANG[i,3]]+y[ANG[i,2]]
    }
  }
  
  # Step 11
  aux=R[ANG[i,2]]
  R[ANG[i,2]]=R[ANG[i,3]]
  R[ANG[i,3]]=aux
  
  iR[R[ANG[i,2]]]=ANG[i,2]
  iR[R[ANG[i,3]]]=ANG[i,3]
  
  
  EXT=rbind(EXT,M[,1:2])
  }
return(EXT)
}

#########################################################################################################
#Checks if a point is an extreme point

check.point <- function(s,x,y,alf,bet,v,iR){
  ind <- FALSE
  X <- cbind(x,y)
  u1 <- c(cos(alf),sin(alf))
  u2 <- c(cos(bet),sin(bet))
  if((sum(X[iR[s],]*u1)<= sum(v[1:2]*u1)) & (sum(v[1:2]*u2) <=sum(X[iR[s+1],]*u2))){ind<-TRUE}
  if((sum(X[iR[s],]*u2)<= sum(v[1:2]*u2)) & (sum(v[1:2]*u1) <=sum(X[iR[s+1],]*u1))){ind<-TRUE}
  return(ind)
}

#########################################################################################################
#Determines the angle of a vector wrt the x-axis
x.angle <- function(u){
  ux <- u[1]
  uy <- u[2]
  angle <- atan2(uy,ux)
  if(ux<0&uy<0)angle <- 2*pi+atan2(uy,ux)
  return(angle)
}

#########################################################################################################
#Determines 1st candidate to be an extreme point

vector_mas <- function(s,x,y,mx,my,alpha,iR){
  n <- length(x)
  ymas <- vector(length=6)
  ymas[5]<-sum(x[iR[(s+1):n]])
  ymas[6]<-sum(y[iR[(s+1):n]])
  ymas[1]<-(1-alpha)/(alpha*n+s*(1-2*alpha))*(n*mx+(2*alpha-1)/(1-alpha)*ymas[5])
  ymas[2]<-(1-alpha)/(alpha*n+s*(1-2*alpha))*(n*my+(2*alpha-1)/(1-alpha)*ymas[6])
  ymas[3]<-s
  ymas[4]<-x.angle(c(ymas[1]-mx,ymas[2]-my))
  return(ymas)
}

#########################################################################################################
#Determines the candidate to be an extreme point when s changes to s+1

vector_masp <- function(s,x,y,mx,my,alpha,iR,y0){
  n <- length(x)
  ymas <- vector(length=6)
  ymas[5]<-y0[5]-x[iR[s]]
  ymas[6]<-y0[6]-y[iR[s]]
  ymas[1]<-(1-alpha)/(alpha*n+s*(1-2*alpha))*(n*mx+(2*alpha-1)/(1-alpha)*ymas[5])
  ymas[2]<-(1-alpha)/(alpha*n+s*(1-2*alpha))*(n*my+(2*alpha-1)/(1-alpha)*ymas[6])
  ymas[3]<-s
  ymas[4]<-x.angle(c(ymas[1]-mx,ymas[2]-my))
  return(ymas)
}

#########################################################################################################
#Determines the candidate to be an extreme point when s changes to s-1

vector_masm <- function(s,x,y,mx,my,alpha,iR,y0){
  n <- length(x)
  ymas <- vector(length=6)
  ymas[5]<-y0[5]+x[iR[s+1]]
  ymas[6]<-y0[6]+y[iR[s+1]]
  ymas[1]<-(1-alpha)/(alpha*n+s*(1-2*alpha))*(n*mx+(2*alpha-1)/(1-alpha)*ymas[5])
  ymas[2]<-(1-alpha)/(alpha*n+s*(1-2*alpha))*(n*my+(2*alpha-1)/(1-alpha)*ymas[6])
  ymas[3]<-s
  ymas[4]<-x.angle(c(ymas[1]-mx,ymas[2]-my))
  return(ymas)
}

#########################################################################################################
#Determines the candidate to be an extreme point for new ordering

vector_mase <- function(s,x,y,mx,my,alpha,iR,y0){
  n <- length(x)
  ymas <- vector(length=6)
  ymas[5]<-y0[5]
  ymas[6]<-y0[6]
  ymas[1]<-(1-alpha)/(alpha*n+s*(1-2*alpha))*(n*mx+(2*alpha-1)/(1-alpha)*ymas[5])
  ymas[2]<-(1-alpha)/(alpha*n+s*(1-2*alpha))*(n*my+(2*alpha-1)/(1-alpha)*ymas[6])
  ymas[3]<-s
  ymas[4]<-x.angle(c(ymas[1]-mx,ymas[2]-my))
  return(ymas)
}

#########################################################################################################
#Determines the angles between pairs of points

xy.angles<-function(x,y){
  # x,y are the coordinates of data points
  n<-length(x)
  m<-choose(n,2)
  M<-matrix(0,m,3)
  
  index <- 1
  dx<-0
  dy<-0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dx<-x[j]-x[i]
      dy<-y[j]-y[i]
      if(dx>=0&dy>=0){
        M[index,]<-c(pi/2+atan2(dy,dx),i,j)}
      else if(dx>=0&dy<0){
        M[index,]<-c(pi/2+atan2(dy,dx),i,j)}
      else if(dx<0&dy>=0){
        M[index,]<-c(pi/2+atan2(-dy,-dx),i,j)}
      else if(dx<0&dy<0){
        M[index,]<-c(pi/2+atan2(-dy,-dx),i,j)}
      index <- index + 1}
  }
  M<-M[order(M[,1]),]
  return(M)
}