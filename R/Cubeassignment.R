
data <- read.csv2("data/cube.csv")
head(data)

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