library(tree)
data <-read.csv2("C:/Users/Dator/Documents/R_HW/ML-Lab-2/data/State.csv")
head(data)

plot(data$MET,data$EX, main=c("per capita State and local expenditures ($)",
                              "vs percentage of people living in metropolitan areas"),
ylab="$ expenditures per capita",
xlab="percentage of people in metropolitan areas",pch="+")

regtree <- tree(data$EX ~ data$MET,data, control = tree.control(48,minsize=2))
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
hist(resid,main=c("Residuals of best regression tree prediction", "of per capita expenditure"),
xlab="residual ($)")
plot(data$MET,pred,col="red",ylim=c(240,400),main="Regression tree fitted values and original values",
     ylab="$ expenditures per capita",
     xlab="percentage of people in metropolitan areas")
points(data$MET,data$EX,pch="+")
legend(x=20,y=400,c("original values","fitted values"),
       pch=c("+","o"),
       col=c("black","red"))

