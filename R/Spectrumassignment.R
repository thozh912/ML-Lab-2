library(XLConnect)
library(fastICA)
library(pls)
#FROM THIS FILE LOCATION, EXCEL FILES SHOULD BE FOUND IN A SUBFOLDER IN THIS FILE LOCATION CALLED DATA
wb = loadWorkbook("D:/R_HW/ML-Lab-2/data/NIRSpectra.xls")
data2 = readWorksheet(wb, sheet = "NIRSpectra", header = TRUE)
data2 <- data2[,-1]
data2 <- data2[,-length(data2)]
data2 <-data2[complete.cases(data2),]
#head(data2)
res=prcomp(data2)
resul <- princomp(data2)
lambda=res$sdev^2
#eigenvalues
# lambda
#proportion of variation
sprintf("%2.3f",lambda/sum(lambda)*100)
# two pcs are needed to explain more than 99% of 
screeplot(res,main="Largest contributions to variance",xlab="Principal components")
plot(res$x[,1],res$x[,2],pch="+",main="PCA scoreplot of PC2 against PC1, unequal axes",xlab="PC1",ylab="PC2",cex=0.5)
# biplot(res)
# U=res$rotation
# head(U)
U=loadings(resul)
plot(U[,1], main="PCA Traceplot, PC1")
plot(U[,2],main="PCA Traceplot, PC2")

res2 <- fastICA(data2,2)
#res2$W is the unmixing matrix "X"W = S where S contain independent components
#res2$K is a pre-whitening matrix which projects data onto first 2 principal components XKW = S
projectmatr <- res2$K %*% res2$W
plot(res2$K[,1], main="ICA K-matrix Traceplot, PC1")
plot(res2$K[,2], main="ICA K-matrix Traceplot, PC2")
#^^These are same as for PCA
plot(projectmatr[,1], main="ICA Traceplot, PC1")
plot(projectmatr[,2],main="ICA Traceplot, PC2")
plot(res2$S,pch="+",main="ICA scoreplot of PC2 against PC1, unequal axes",xlab="PC1",ylab="PC2",cex=0.5)


data2 <- readWorksheet(wb, sheet = "NIRSpectra", header = TRUE)
data2 <- data2[,-1]
data2 <-data2[complete.cases(data2),]
n <- dim(data2)[1]
set.seed(12345)
ind <- sample(1:n,floor(0.5*n))
train <- data2[ind,]
test <- data2[-ind,]

set.seed(12345)
pcr.fit=pcr(Viscosity~., data=train, validation="CV") 
validationplot(pcr.fit,val.type="MSEP",main="MSEP scores for CV:ed PCR models") 
pcrmodel=pcr(Viscosity~., 20,data=train, validation="none") 
pcr.pred <- predict(pcrmodel, newdata=test, ncomp = 20) 
pcr.mse <- 1/n * sum((test$Viscosity - pcr.pred)^2, na.rm=TRUE)
paste("MSE for test data for PCR model with 20 components:", signif(pcr.mse,3)) 

set.seed(12345)
plsr.fit <- plsr(Viscosity~., data=train, validation="CV")
validationplot(plsr.fit,val.type="MSEP",main="MSEP scores for CV:ed PLS models") 
plsmodel=plsr(Viscosity~., 10,data=train, validation="none") 
pls.pred <- predict(plsmodel, newdata=test, ncomp = 10) 
pls.mse <- 1/n * sum((test$Viscosity - pls.pred)^2, na.rm=TRUE) 
paste("MSE for test data for PLS model with 10 components:", signif(pls.mse,3)) 
