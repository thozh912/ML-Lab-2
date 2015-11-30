
data <- read.csv2("D:/R_HW/ML-Lab-2/data/cube.csv")
#head(data)

piecew_const <- function(resp,input,knots){
  intervals <- list(c(min(input)-0.01,knots[1]))
  if(length(knots) > 1){
    for(i in 1:(length(knots)-1)){
      intervals[[i+1]] <- c(knots[i],knots[i+1])
    }
    intervals[[length(knots)+1]] <- c(knots[length(knots)],max(input)+0.01)
  } else{
    intervals[[2]] <- c(knots[1],max(input))
  }

#   sortedx <- sort(data$input,index.return = TRUE)
#   sortedy <- data$resp[sortedx$ix]
  listofgroups <- list()
  
  
  for(i in 1:length(intervals)){
     
    for(j in 1:length(input)){
      if(input[j] > intervals[[i]][1] && input[j] <
         intervals[[i]][2]){
        if(length(listofgroups) < i){
          listofgroups[[i]] <- resp[j]
        } else{
        listofgroups[[i]] <- c(listofgroups[[i]],resp[j])
        }
      }
    }
  }
  
  means <- c()
  for(i in 1:length(listofgroups)){
    means[i] <- mean(listofgroups[[i]])
  }
  plot(input,resp)
  for(i in 1:length(knots)){
    abline(v=knots[i],lty=2,col="gray")
    lines(intervals[[i]],rep(means[i],2),col="red")
  }
  lines(intervals[[length(knots)+1]],rep(means[length(knots)+1],2),col="red")
  return(means)
}

means <- piecew_const(data$y,data$x,c(3,6))
library(mgcv)
data <- read.csv2("D:/R_HW/ML-Lab-2/data/influenza.csv")

influ_data_ts <- ts(data$Influenza,start=c(data$Year[1],
                                     data$Week[1]),freq=52)
mort_data_ts <- ts(data$Mortality,start=c(data$Year[1],
                                           data$Week[1]),freq=52)
plot(influ_data_ts,main="Weekly Influenza cases in Sweden vs year",
     xlab="year",ylab="Confirmed Influenza cases")
plot(mort_data_ts,main="Weekly Mortality in Sweden vs year",
     xlab="year",ylab="Mortality")
acf(influ_data_ts,
    main="Auto-correlation function of weekly influenza cases",
    xlab="Years lag",
    lag.max=52*3)
acf(mort_data_ts,
    main="Auto-correlation function of weekly Mortality",
    xlab="Years lag",
    lag.max=52*3)


gam_mort <- gam(Mortality ~ Year + s(Week),data=data)

#str(gam_mort)
plot(data$Mortality,type="l",
     main=c("Weekly Mortality in Sweden vs year with GAM fit",
            "linear component for year,spline function for week"),
     xlab="# weeks after week 52 1994",ylab="Mortality")
lines(gam_mort$fitted.values,col="red")
summary(gam_mort)

plot(gam_mort, main="Week spline component over a year")
gam_mort2 <- gam(formula=Mortality ~ Year + s(Week, k=3), data=data) 
gam_mort3 <- gam(formula=Mortality ~ Year + s(Week, k=20), data=data) 



plot(data$Mortality,type="l",
     main=c("Weekly Mortality in Sweden vs year with GAM fit",
     "week spline function space has dimension 3"),
     xlab="# weeks after week 52 1994",ylab="Mortality")
lines(gam_mort2$fitted.values, col="red") 
paste("Deviance explained by fit, Week spline basis dimension 3:",61.6,"%")
plot(data$Mortality,type="l",
     main=c("Weekly Mortality in Sweden vs year with GAM fit",
            "week spline function space has dimension 20"),
     xlab="# weeks after week 52 1994",ylab="Mortality")
lines(gam_mort3$fitted.values, col="red") 
paste("Deviance explained by fit, Week spline basis dimension 20:",68.8,"%")
plot(data$Time,data$Influenza,ylim=c(-250,600),
     main="Weekly influenza and residuals of mortality GAM fit",
     xlab="year")
points(data$Time,gam_mort$residuals,col="red")
legend("topleft",legend=c("influenza cases","Residual of GAM fit"),
       pch=c("o","o"),col=c("black","red"))

gam_mort4 <- gam(formula=Mortality ~ s(Week, k=5) + s(Influenza, k=4), data=data)
summary(gam_mort4)
par(mfrow=c(1,2)) 
plot(gam_mort4) 
par(mfrow=c(1,1)) 



plot(data$Mortality,type="l",
     main=c("Weekly Mortality in Sweden vs year with GAM fit",
            "spline functions of week, year and influenza cases"),
     xlab="# weeks after week 52 1994",ylab="Mortality")
lines(gam_mort4$fitted.values, col="red") 
paste("Deviance explained by this fit:",70.9,"%")
sse <- sum((data$Mortality-gam_mort4$fitted.values)^2)
paste("SSE for best GAM fit:",round(sse,0))
## NA
