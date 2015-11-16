library(tree)
library(boot)
data <-read.csv2("C:/Users/Dator/Documents/R_HW/ML-Lab-2/data/State.csv")
data <- data[order(data$MET),]
plot(data$MET,data$EX, main=c("per capita State and local expenditures ($)",
                              "vs percentage of people living in metropolitan areas"),
ylab="$ expenditures per capita",
xlab="percentage of people in metropolitan areas",pch="+")

regtree <- tree(EX ~ MET,data, control = tree.control(48,minsize=2))
set.seed(12345)
cvresult <- cv.tree(regtree)
# best result at size = 3
plot(cvresult$size,cvresult$dev, type="b", main="Deviance of fitted tree against tree size",
     xlab="Size",ylab="Deviance")

bestregtree <- prune.tree(regtree,best=3)
plot(bestregtree)
text(bestregtree,pretty=1)

pred <- predict(bestregtree,data)

resid <- data$EX - pred
#residuals look uniform

hist(resid,main=c("Residuals of best regression tree prediction", "of per capita expenditure"),
xlab="residual ($)")
plot(data$MET,pred,col="red",ylim=c(240,400),main="Regression tree fitted values and original values",
     ylab="$ expenditures per capita",
     xlab="percentage of people in metropolitan areas")
points(data$MET,data$EX,pch="+")
legend(x=20,y=400,c("original values","fitted values"),
       pch=c("+","o"),
       col=c("black","red"))
f <- function(datainp,ind){
  data1 <- datainp[ind,]
  res <- tree(EX ~ MET,data1, control = tree.control(dim(data1)[1],minsize=2))
  bestrestree <- prune.tree(res,best=3)
  predictions <- predict(bestrestree, newdata=data)
  return(predictions)
}

set.seed(12345)
bootobj1 <- boot(data,f, R=1000)
confintvs <- envelope(bootobj1)
#plot(bootobj1)

plot(data$MET,pred,col="red",ylim=c(150,500),
     main=c("Regression tree fitted values and original values",
            "and 95% non-parametric bootstrap confidence bands"),
     ylab="$ expenditures per capita",
     xlab="percentage of people in metropolitan areas")
points(data$MET,data$EX,pch="+")
points(data$MET,confintvs$point[2,], type="l", col="blue")
points(data$MET,confintvs$point[1,], type="l", col="blue")

legend(x=20,y=500,c("original values","fitted values","confidence intervals"),
       pch=c("+","o",NA),lwd=1,lty=c(NA,NA,1),
       col=c("black","red","blue"))
rng <- function(data2,mle){
  data1 = data.frame(MET = data2$MET, EX = data2$EX)
  n = length(data2$EX)
  data1$EX = rnorm(n,predict(mle, newdata=data1),
                   sd(data$EX-predict(mle, newdata=data1)))
  return(data1)
}

f1 = function(data1){
  res <- tree(EX ~ MET,data1, control = tree.control(dim(data1)[1],minsize=2))
  bestrestree <- prune.tree(res,best=3)
  expenditures <- predict(bestrestree, newdata=data)
  return(expenditures)
}

f2 = function(data1){
  res <- tree(EX ~ MET,data1, control = tree.control(dim(data1)[1],minsize=2))
  bestrestree <- prune.tree(res,best=3)
  expenditures <- rnorm(dim(data)[1],predict(bestrestree, newdata=data),
                        sd(resid))
  return(expenditures)
}

set.seed(12345)
bootobj2 <- boot(data, statistic= f1,R=1000,
                 mle=bestregtree,ran.gen=rng,sim="parametric")
set.seed(12345)
bootobj3 <- boot(data,statistic= f2,R=1000,
                 mle=bestregtree,ran.gen=rng,sim="parametric")
#plot(bootobj2)
confintvs2 <- envelope(bootobj2)
confintvs3 <- envelope(bootobj3)
plot(data$MET,pred,col="red",ylim=c(150,550),
     main=c("Regression tree fitted values and original values",
            "95% parametric bootstrap confidence and prediction bands"),
     ylab="$ expenditures per capita",
     xlab="percentage of people in metropolitan areas")
points(data$MET,data$EX,pch="+")
points(data$MET,confintvs2$point[2,], type="l", col="blue")
points(data$MET,confintvs2$point[1,], type="l", col="blue")
points(data$MET,confintvs3$point[2,], type="l", col="pink")
points(data$MET,confintvs3$point[1,], type="l", col="pink")
legend(x=20,y=550,c("original values","fitted values",
                    "confidence bands","prediction bands"),
       pch=c("+","o",NA,NA),lwd=1,lty=c(NA,NA,1,1),
       col=c("black","red","blue","pink"))

library(XLConnect)
library(fastICA)
library(pls)
#FROM THIS FILE LOCATION, EXCEL FILES SHOULD BE FOUND IN A SUBFOLDER IN THIS FILE LOCATION CALLED DATA
wb = loadWorkbook("C:/Users/Dator/Documents/R_HW/ML-Lab-2/data/NIRSpectra.xls")
data2 = readWorksheet(wb, sheet = "NIRSpectra", header = TRUE)
data2 <- data2[,-1]
data2 <- data2[,-length(data2)]
#data2 <-data2[complete.cases(data2),]
#head(data2)
set.seed(12345)
res=prcomp(data2)
#resul <- princomp(data2,scale=TRUE)
lambda=res$sdev^2
percentage <- lambda/sum(lambda)*100
#eigenvalues
# lambda
#proportion of variation
paste("Percentage of variance explained by first PC: ",signif(percentage[1],3))
paste("Percentage of variance explained by second PC: ",signif(percentage[2],2))
# two pcs are needed to explain 99% of 
screeplot(res,main="Largest contributions to variance",xlab="Principal components")
plot(res$x[,1],res$x[,2],pch="+",
     main="PCA scoreplot of PC2 against PC1, unequal axes",
     xlab="PC1",ylab="PC2",cex=0.5)
# biplot(res)
# U=res$rotation
# head(U)
#U=loadings(resul)
plot(res$rotation[-nrow(res$rotation),1],
     main="PCA Traceplot, PC1",
     xlab="component",ylab="coefficients")
plot(res$rotation[-nrow(res$rotation),2],
     main="PCA Traceplot, PC2",
     xlab="component",ylab="coefficients")
# plot(U[,1], main="PCA Traceplot, PC1")
# plot(U[,2],main="PCA Traceplot, PC2")
set.seed(12345)
res2 <- fastICA(data2,2)
#res2$W is the unmixing matrix "X"W = S where S contain independent components
#res2$K is a pre-whitening matrix which projects data onto first 2 principal components XKW = S
projectmatr <- res2$K %*% res2$W
#plot(res2$K[,1], main="ICA K-matrix Traceplot, PC1")
#plot(res2$K[,2], main="ICA K-matrix Traceplot, PC2")
#^^These are same as for PCA
plot(projectmatr[,1], main="ICA Traceplot, PC1")
plot(projectmatr[,2],main="ICA Traceplot, PC2")
plot(res2$S,pch="+",
     main="ICA scoreplot of PC2 against PC1, unequal axes",
     xlab="PC1",ylab="PC2",cex=0.5)

data2 <- readWorksheet(wb, sheet = "NIRSpectra", header = TRUE)
data2 <- data2[,-1]
#data2 <-data2[complete.cases(data2),]
set.seed(12345)
ind <- sample(1:784,392)
train <- data2[ind,]
test <- data2[-ind,]
n <- dim(test[complete.cases(test),])[1]

set.seed(12345)
pcr.fit=pcr(Viscosity~., data=train, validation="CV") 
validationplot(pcr.fit,val.type="MSEP",main="MSEP scores for CV:ed PCR models") 
pcrmodel=pcr(Viscosity~., 20,data=train, validation="none") 
pcr.pred <- predict(pcrmodel, newdata=test[complete.cases(test),], ncomp = 20) 
pcr.mse <- 1/n * sum((test[complete.cases(test),]$Viscosity - pcr.pred)^2, na.rm=TRUE)
paste("MSE for test data for PCR model with 20 components:", signif(pcr.mse,3)) 

set.seed(12345)
plsr.fit <- plsr(Viscosity~., data=train, validation="CV")
validationplot(plsr.fit,val.type="MSEP",main="MSEP scores for CV:ed PLS models") 
plsmodel=plsr(Viscosity~., 10,data=train, validation="none") 
pls.pred <- predict(plsmodel, newdata=test[complete.cases(test),], ncomp = 10) 
pls.mse <- 1/n * sum((test[complete.cases(test),]$Viscosity - pls.pred)^2, na.rm=TRUE) 
paste("MSE for test data for PLS model with 10 components:", signif(pls.mse,3)) 
## NA
