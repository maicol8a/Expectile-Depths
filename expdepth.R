expdepth <- function(x,data){
  if(min(chull(rbind(x,data)))==1) {return(0)}
  n <- nrow(data)
  # Step 1
  data <- t(t(data)-x)
  # Step 2
  sum.data <- colSums(data)
  # Step 3
  gamma <- atan2(data[,2],data[,1])
  # Steps 4,5
  ANG <- cbind(data,gamma+pi*((gamma<(-pi/2))-(gamma>pi/2)),
               2*(gamma>=-pi/2)*(gamma<=pi/2)-1)
  ord <- order(ANG[,3])
  ANG <- ANG[ord,]
  depth <- 1e16
  # Steps 6
  sum.pos <- colSums(ANG[,1:2]*(ANG[,4]==1))
  sum.neg <- colSums(ANG[,1:2]*(ANG[,4]==-1))
  for(i in 1:n) {
    sum.pos <- sum.pos+ANG[i,1:2]*(ANG[i,4]==-1)
    sum.neg <- sum.neg+ANG[i,1:2]*(ANG[i,4]==1)
    # Step 8
    depth1 <- (sum.data[1]*cos(ANG[i,3]+pi/2)+sum.data[2]*sin(ANG[i,3]+pi/2))/(sum.pos[1]*cos(ANG[i,3]+pi/2)+sum.pos[2]*sin(ANG[i,3]+pi/2))
    depth2 <- (sum.data[1]*cos(ANG[i,3]-pi/2)+sum.data[2]*sin(ANG[i,3]-pi/2))/(sum.neg[1]*cos(ANG[i,3]-pi/2)+sum.neg[2]*sin(ANG[i,3]-pi/2))
    depth <- min(depth,depth1,depth2)
    # Step 9
    sum.pos <- sum.pos-ANG[i,1:2]*(ANG[i,4]==1)
    sum.neg <- sum.neg-ANG[i,1:2]*(ANG[i,4]==-1)
  }
  # Step 10
  expdepth <- 1/(2-depth)
  return(expdepth)
}
