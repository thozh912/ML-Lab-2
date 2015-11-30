library(MASS)

www = "C:/Users/Dator/Documents/R_HW/ML-Lab-2/data/JapanTemp.dat"
data <- read.delim(www, header = TRUE, sep="")

Kernel <- function(x_0,x,lambda,sigmaf){
  res <- matrix(0,nrow=length(x_0),ncol=length(x))
  for(i in 1:length(x_0)){
    for(j in 1:length(x)){
      res[i,j] <- sigmaf^2 * exp(-1/2 * (abs(x_0[i] - x[j]) / lambda )^2)
    }
  }
  return(res)
}

posteriorGP <- function(x,y,xStar,hyperParam,sigmaNoise){
  
  meanfstar <- Kernel(xStar,x,hyperParam[2],hyperParam[1]) %*% 
    solve(Kernel(x,x,hyperParam[2],hyperParam[1]) +
    sigmaNoise^2 * diag(length(x))) %*% y
  covfstar <- Kernel(xStar,xStar,hyperParam[2],hyperParam[1]) -
    Kernel(xStar,x,hyperParam[2],hyperParam[1]) %*% 
    solve(Kernel(x,x,hyperParam[2],hyperParam[1]) +
    sigmaNoise^2 * diag(length(x))) %*%
    Kernel(x,xStar,hyperParam[2],hyperParam[1])
  posterior <- mvrnorm(1,meanfstar,covfstar)
  
  return(list(posteriorMean = meanfstar, posteriorCov= covfstar))
}

btask <- posteriorGP(0.4,0.719,seq(-1,1,0.01),c(1, 0.3),0.1)
uband <- btask$posteriorMean + 2 * sqrt(diag(btask$posteriorCov))
lband <- btask$posteriorMean - 2 * sqrt(diag(btask$posteriorCov))

plot(seq(-1,1,0.01),ylim=c(-4,4),btask$posteriorMean,type="l",
     main=c("Gaussian Process Regression, one observation",
            "with 95% pointwise confidence intervals"),ylab="mean(posterior distribution)",
     xlab="x")
legend("top",c(expression(paste(lambda," = 0.3")),
               expression(paste(sigma[f]," = 1")),
               expression(paste(sigma[n]," = 0.1"))),bty="n")
     
lines(seq(-1,1,0.01),uband,col="red",lty=2)
lines(seq(-1,1,0.01),lband,col="red",lty=2)
ctask <- posteriorGP(c(0.4,-0.6),c(0.719,0.044),seq(-1,1,0.01),c(1, 0.3),0.1)
uband2 <- ctask$posteriorMean + 2 * sqrt(diag(ctask$posteriorCov))
lband2 <- ctask$posteriorMean - 2 * sqrt(diag(ctask$posteriorCov))
plot(seq(-1,1,0.01),ylim=c(-3,3),ctask$posteriorMean,type="l",
     main=c("Gaussian Process Regression, two observations",
            "with 95% pointwise confidence intervals"),ylab="mean(posterior distribution)",
      xlab="x")
legend("top",c(expression(paste(lambda," = 0.3")),
               expression(paste(sigma[f]," = 1")),
               expression(paste(sigma[n]," = 0.1"))),bty="n")
lines(seq(-1,1,0.01),uband2,col="red",lty=2)
lines(seq(-1,1,0.01),lband2,col="red",lty=2)
yvec <- c(0.768,-0.044,-0.940,0.719,-0.664)
xvec <- c(-1,-.6,-.2,.4,.8)
dtask <- posteriorGP(xvec,yvec,seq(-1,1,0.01),c(1, 0.3),0.1)
uband3 <- dtask$posteriorMean + 2 * sqrt(diag(dtask$posteriorCov))
lband3 <- dtask$posteriorMean - 2 * sqrt(diag(dtask$posteriorCov))
plot(seq(-1,1,0.01),ylim=c(-3,3),dtask$posteriorMean,type="l",
     main=c("Gaussian Process Regression, five observations",
            "with 95% pointwise confidence intervals"),ylab="mean(posterior distribution)",
     xlab="x")
legend("top",c(expression(paste(lambda," = 0.3")),
               expression(paste(sigma[f]," = 1")),
               expression(paste(sigma[n]," = 0.1"))),bty="n")
lines(seq(-1,1,0.01),uband3,col="red",lty=2)
lines(seq(-1,1,0.01),lband3,col="red",lty=2)

etask <- posteriorGP(xvec,yvec,seq(-1,1,0.01),c(1, 1),0.1)
uband4 <- etask$posteriorMean + 2 * sqrt(diag(etask$posteriorCov))
lband4 <- etask$posteriorMean - 2 * sqrt(diag(etask$posteriorCov))
plot(seq(-1,1,0.01),ylim=c(-3,3),etask$posteriorMean,type="l",
     main=c("Gaussian Process Regression, five observations",
            "with 95% pointwise confidence intervals"),ylab="mean(posterior distribution)",
     xlab="x")
legend("top",c(expression(paste(lambda," = 1")),
               expression(paste(sigma[f]," = 1")),
               expression(paste(sigma[n]," = 0.1"))),bty="n")
lines(seq(-1,1,0.01),uband4,col="red",lty=2)
lines(seq(-1,1,0.01),lband4,col="red",lty=2)

jtask <- posteriorGP(data$time,data$temp,seq(0,1,0.01),c(1, 1),0.1)
ubandj <- jtask$posteriorMean + 2 * sqrt(diag(jtask$posteriorCov))
lbandj <- jtask$posteriorMean - 2 * sqrt(diag(jtask$posteriorCov))

plot(seq(0,1,0.01),jtask$posteriorMean,type="l", ylim = c(10,35),
     main=c("Gaussian Process Regression","Japan temperature over a year"
            ),ylab="mean temp (posterior distribution)",
     xlab="time fraction of full year")
points(data$time,data$temp,pch="+",col="gray",cex=0.6)
legend("topleft",c(expression(paste(lambda," = 1")),
               expression(paste(sigma[f]," = 1")),
               expression(paste(sigma[n]," = 0.1"))),bty="n")


jtask <- posteriorGP(data$time,data$temp,seq(0,1,0.01),c(1, 0.3),0.1)
ubandj <- jtask$posteriorMean + 2 * sqrt(diag(jtask$posteriorCov))
lbandj <- jtask$posteriorMean - 2 * sqrt(diag(jtask$posteriorCov))

plot(seq(0,1,0.01),jtask$posteriorMean,type="l", ylim = c(10,35),
     main=c("Gaussian Process Regression","Japan temperature over a year"
     ),ylab="mean temp (posterior distribution)",
     xlab="time fraction of full year")
points(data$time,data$temp,pch="+",col="gray",cex=0.6)
legend("topleft",c(expression(paste(lambda," = 0.3")),
                   expression(paste(sigma[f]," = 1")),
                   expression(paste(sigma[n]," = 0.1"))),bty="n")

jtask <- posteriorGP(data$time,data$temp,seq(0,1,0.01),c(1, 1),1)
ubandj <- jtask$posteriorMean + 2 * sqrt(diag(jtask$posteriorCov))
lbandj <- jtask$posteriorMean - 2 * sqrt(diag(jtask$posteriorCov))

plot(seq(0,1,0.01),jtask$posteriorMean,type="l", ylim = c(10,35),
     main=c("Gaussian Process Regression","Japan temperature over a year"
     ),ylab="mean temp (posterior distribution)",
     xlab="time fraction of full year")
points(data$time,data$temp,pch="+",col="gray",cex=0.6)
legend("topleft",c(expression(paste(lambda," = 1")),
                   expression(paste(sigma[f]," = 1")),
                   expression(paste(sigma[n]," = 1"))),bty="n")

jtask <- posteriorGP(data$time,data$temp,seq(0,1,0.01),c(0.05, 1),1)

plot(seq(0,1,0.01),jtask$posteriorMean,type="l", ylim = c(10,35),
     main=c("Gaussian Process Regression","Japan temperature over a year"
     ),ylab="mean temp (posterior distribution)",
     xlab="time fraction of full year")
points(data$time,data$temp,pch="+",col="gray",cex=0.6)
legend("topleft",c(expression(paste(lambda," = 1")),
                   expression(paste(sigma[f]," = 0.05")),
                   expression(paste(sigma[n]," = 1"))),bty="n")

## NA
