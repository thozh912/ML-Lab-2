library(mgcv)
data <- read.csv2("data/influenza.csv")

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
summary(gam_mort)
#str(gam_mort)
plot(data$Mortality,type="l",
     main=c("Weekly Mortality in Sweden vs year with GAM fit",
            "linear component for year,spline function for week"),
     xlab="# weeks after week 52 1994",ylab="Mortality")
lines(gam_mort$fitted.values,col="red")

plot(gam_mort, main="Week spline component over a year")

gam_mort2 <- gam(formula=Mortality ~ Year + s(Week, k=3), data=data) 
gam_mort3 <- gam(formula=Mortality ~ Year + s(Week, k=20), data=data) 



plot(data$Mortality,type="l",
     main=c("Weekly Mortality in Sweden vs year with GAM fit",
     "week spline function space has dimension 3"),
     xlab="# weeks after week 52 1994",ylab="Mortality")
lines(gam_mort2$fitted.values, col="red") 
paste("Deviance of Mortality fit, Week spline basis dimension 3:",round(gam_mort2$deviance,0))
plot(data$Mortality,type="l",
     main=c("Weekly Mortality in Sweden vs year with GAM fit",
            "week spline function space has dimension 20"),
     xlab="# weeks after week 52 1994",ylab="Mortality")
lines(gam_mort3$fitted.values, col="red") 
paste("Deviance of Mortality fit, Week spline basis dimension 20:",round(gam_mort3$deviance,0))

plot(data$Time,data$Influenza,ylim=c(-250,400),
     main="Weekly influenza and residuals of mortality GAM fit",
     xlab="year")
points(data$Time,gam_mort$residuals,col="red")
legend("topleft",legend=c("influenza cases","Residual of GAM fit"),
       pch=c("o","o"),col=c("black","red"))

gam_mort4 <- gam(formula=Mortality ~ s(Week, k=5) + s(Influenza, k=4), data=data)

par(mfrow=c(1,2)) 
plot(gam_mort4) 
par(mfrow=c(1,1)) 

summary(gam_mort4)

plot(data$Mortality,type="l",
     main=c("Weekly Mortality in Sweden vs year with GAM fit",
            "spline functions of week, year and influenza cases"),
     xlab="# weeks after week 52 1994",ylab="Mortality")
lines(gam_mort4$fitted.values, col="red") 

sse <- sum((data$Mortality-gam_mort4$fitted.values)^2)
paste("SSE for best GAM fit:",round(sse,0))
