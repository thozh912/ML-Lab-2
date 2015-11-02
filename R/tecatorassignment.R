library(XLConnect)
library(glmnet)
library(MASS)
#hehehehehe
#FROM THIS FILE LOCATION, EXCEL FILES SHOULD BE FOUND IN A SUBFOLDER IN THIS FILE LOCATION CALLED DATA
wb = loadWorkbook("data/tecator.xlsx")
data2 = readWorksheet(wb, sheet = "data", header = TRUE)

plot(data2$Protein,data2$Moisture,pch="+",main="Plot of Moisture against Protein",
     xlab = "Protein", ylab = "Moisture")

# $y_{j} = \sum_{i=0}^{n} w_{i}x_{j}^{i} + e$ , $e$ is $N(0,\sigma^2)$ distributed.

n=dim(data2)[1]
set.seed(12345)
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

#AIC = -2loglik(D) + 2d(model) where d(model) is the no. of parameters. sigma2 counts as a parameter
# We make assumption that the outcomes are normally distributed where for observation j mean is 
# \sum_{i=0}^{n} w_{i}x_{j}^{i} and variance is \sigma^2

#We estimate \sigma^2 by the MSE of linear model of total data set.

Xlinearforentire <- cbind(rep(1,n),data2$Protein)
predlinearforentire <- Xlinearforentire %*% wvalues[[1]]
sigma2 <- sum(((data2$Moisture) - predlinearforentire)^2) / n

AIC_calc <- function(data2,wvalues,sigma2){
  AICs <- c()
  
  for(j in 1:6){
    jpt <- 1
    Xentire <- rep(1,n)
    while(jpt <= j){
      Xentire <- cbind(Xentire,data2$Protein^jpt)
      jpt <- jpt + 1
    }
    entirepredicted <- Xentire %*% wvalues[[j]]
    loglik <- n * (log(2 * pi) + log(sigma2)) + 
      1 / sigma2 * sum((data2$Moisture - entirepredicted)^2)
    AIC <- loglik + 2 * (j + 2)
    AICs <- c(AICs , AIC)
  }
  return(AICs)
}

AICvec <- AIC_calc(data2,wvalues,sigma2)
for(j in 1:6){
  print(paste("AIC at polynomial degree",j,"is :",AICvec[j]))
}

fat.lm <- lm(Fat ~.,data = data2[,2:102])
fat.lm2 <- stepAIC(fat.lm,trace = FALSE)
fat.lm2$anova

covariates <- scale(data2[,2:101])
response <- scale(data2[,102])
ridgemod <- glmnet(as.matrix(covariates), response,
                   family = "gaussian", alpha = 0)
plot(ridgemod, xvar = "lambda" , label = TRUE)

lassomod <- glmnet(as.matrix(covariates), response,
                   family = "gaussian", alpha = 1)
plot(lassomod, xvar = "lambda" , label = TRUE)

#lambda search for smaller lambda is not very fruitful since the sd of cvm
#is larger than the improvement.
lassocv <- cv.glmnet(as.matrix(covariates), response, alpha =1)
plot(lassocv)
lassocv$lambda.min
coef(lassocv, s= "lambda.min")
# 14 parameters are non-zero, including inercept