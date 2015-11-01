library(XLConnect)
library(glmnet)

#heheheheheke
#FROM THIS FILE LOCATION, EXCEL FILES SHOULD BE FOUND IN A SUBFOLDER IN THIS FILE LOCATION CALLED DATA
wb = loadWorkbook("data/tecator.xlsx")
data2 = readWorksheet(wb, sheet = "data", header = TRUE)

plot(data2$Protein,data2$Moisture,pch="+",main="Plot of Moisture against Protein",
     xlab = "Protein", ylab = "Moisture")

# $y_{j} = \sum_{i=0}^{n} w_{i}x_{j}^{i} + e$ , $e$ is $N(0,\sigma^2)$ distributed.

n=dim(data2)[1]
set.seed(78)
id=sample(1:n, floor(n*0.5))
train=data2[id,]
valid=data2[-id,]

powerlist <-list()
validlist <- list()
wvalues <- list()
for(i in 1:6){
  powerlist[[i]] <- train$Protein^i
  validlist[[i]] <- valid$Protein^i
}

MSEtrain <-c()
MSEvalid <-c()


for(j in 1:6){
  jpt <- 1
  X <- rep(1,dim(train)[1])
  Xvalid <- rep(1,dim(valid)[1])
  while(jpt <= j){
  X <- cbind(X,powerlist[[jpt]])
  Xvalid <- cbind(Xvalid,validlist[[jpt]])
  jpt <- jpt + 1
  }
  
  #Using SVD decomp to find least squares as the conditionnumber is too high 
  
  svdobj <- svd(X)
  #d <- diag(svdobj$d)
  z <- t(svdobj$u) %*% (train$Moisture)
  z <- z / svdobj$d
  wvalues[[j]] <- svdobj$v %*% z
  
 
  #wvalues <- solve(t(X) %*% X) %*% t(X) %*% (train$Moisture)
  
  trainpredicted <- X %*% wvalues[[j]]
  validpredicted <- Xvalid %*% wvalues[[j]]
  
  MSEtrain <- c(MSEtrain,sum(((train$Moisture) - trainpredicted)^2) / dim(train)[1])
  MSEvalid <- c(MSEvalid,sum(((valid$Moisture) - validpredicted)^2) / dim(valid)[1])
}

plot(MSEtrain,type="l",main="MSE for polynomial fit",xlab="Polynomial degree",
     ylab="MSE",ylim=c(20,50),lwd=2.5)
lines(MSEvalid,col="red",lwd=2.5)
legend("topright",c("training data","validation data"),
       lty=c(1,1),
       lwd=c(2.5,2.5),col=c("black","red"))

#AIC = -2loglik(D) + 2d(model) where d(model) is the no. of parameters.
# We make assumption that the outcomes are normally distributed where for observation j mean is 
# \sum_{i=0}^{n} w_{i}x_{j}^{i} and variance is \sigma^2

AIC_calc <- function()