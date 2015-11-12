library(MASS)

data <- longley

# I do not like picking Nfolds = 10 when nrow(X) = 16

Ridgeregression <- function(X,Y,Lambda,Nfolds){
  
  n <- Nfolds
  rows <- 1:nrow(X)
 
  partsize <- floor( nrow(X) / n)
  
  test_index <-list()
  
  set.seed(-56789)
  for(i in 1:n){
    test_index[[i]] <- sample(rows,partsize)
    
    rows <- setdiff(rows,test_index[[i]])
  }
  
  for(k in 1:length(rows)){
    test_index[[k]] <- c(test_index[[k]], rows[k])
  }
  
  cvscore <- 0
  
  for(j in 1:n){
    Xjth <- as.matrix(X[-test_index[[j]],])
    Xjthcenter <- colMeans(Xjth)
    Xjth <- Xjth - matrix(rep(Xjthcenter,nrow(Xjth)),ncol = ncol(Xjth), byrow=TRUE)
    
    Yjth <- as.matrix(Y[-test_index[[j]]])
    Yjthcenter <- colMeans(Yjth)
    Yjth <- Yjth - Yjthcenter
    w_ridge <- solve(t(Xjth) %*% Xjth + Lambda * diag(ncol(Xjth))) %*% t(Xjth) %*% Yjth
    
    jth_pred <- (as.matrix(X[test_index[[j]],]) -
                   matrix(rep(Xjthcenter,nrow(X[test_index[[j]],])),ncol = ncol(Xjth), byrow=TRUE))  %*% w_ridge
   
    loss <- (Y[test_index[[j]]] - Yjthcenter - jth_pred)^2 
    
    loss_sum <- sum(loss)
    
    cvscore <- cvscore + loss_sum 
  }
  
  return(cvscore)
}



scores <- c()
for(lambda in seq(from=1, to=7, by=1)){
  scores <- c(scores,Ridgeregression(data[,1:6],data$Employed,lambda,10))
  
}
which.min(scores)
plot(scores,type="l",ylab="CV-scores",xlab=expression(lambda),
     main=c("CV-score vs",expression(lambda), "value"))
