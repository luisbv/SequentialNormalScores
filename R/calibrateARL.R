source("ARL.r")
n = 5 #subgroup size
m = 100 #reference-sample size
dist = "Normal"
mu = c(0, 0) # c(reference sample mean, monitoring sample mean)
sigma = c(1, 1) # c(reference sample sd, monitoring sample sd)

#### Distribution parameters
dist.par = c(0, 1, 1) # c(location, scale, shape)

#### Control chart parameters
#chart.par = c(3); chart = "Shewhart"
#Shewhart c(k), where k comes from UCL = mu + k*sigma, LCL = mu - k*sigma.

#chart.par = c(0.5, 4, 1); chart = "CUSUM"
#CUSUM c(k, h, t), where k threshold and h is the decision limit. t is the type of the chart (1:positive, 2:negative, 3:two sides)

chart.par = c(0.2, 2.962); chart = "EWMA"
#EWMA c(lambda, L), where lambda is the smoothing constant and L multiplies standard deviation to get control limit

#### Other Parameters
replicates = 500

calibrateARL <- function(targetARL=NULL, targetmRL=NULL, n, m, theta=NULL, Ftheta=NULL, dist, mu, sigma, dist.par=c(0, 1, 1), initial.par, replicates=50000, chart){
  #Check for errors
  if(is.null(targetARL) && is.null(targetmRL)){
    print("ERROR: Target ARL or target mRL missing")
    return()
  }else if(!is.null(targetARL) && !is.null(targetmRL)){
    print("ERROR: Two targets defined, delete one")
    return()
  }
  
  
  switch(chart,
    Shewhart = {
      name.par = "k"
      index.par = 1
    },
    CUSUM = {
      name.par = "h"
      index.par = 2
    },
    EWMA = {
      name.par = "L"
      index.par = 2
    }
  )

  x = rep(NA, 3)
  y = x
  
  x[1] = initial.par[index.par]
  for(i in 1:3){
    initial.par[index.par] = x[i]
    result = getARLSNS(n, m, theta=NULL, Ftheta=NULL, dist, mu, sigma, dist.par=dist.par, chart.par=initial.par, replicates=replicates, chart=chart)
    
    if(!is.null(targetARL)){
      y[i] = result$ARL
      target = targetARL
      name = "ARL"
    }else{
      y[i] = result$mRL
      target = targetmRL
      name = "mRL"
    }
    
    if(abs(y[i] - target) <= 0.05*target){
      cat("Convergence found with", name.par, "=",x[i],"--",name, "=", y[i], "\n", sep=" ")
      output = list(
        objective.function = y[i],
        par.value = x[i]
      )
      return(output)
    }else{
      if(i == 2){
        x0 = x[i-1]
        x1 = x[i]
        y0 = y[i-1]
        y1 = y[i]
        
        m = (y1-y0)/(x1-x0)
        b = y0-m*x0
        x2 = (target-b)/m
        bounded = FALSE
        
        if((x2 - max(c(x1, x0))) > abs(x1 - x0)){
          x2 = max(c(x1, x0)) + abs(x1 - x0)
          bounded = TRUE
        }else if((min(c(x1, x0)) - x2) > abs(x1 - x0)){
          x2 = min(c(x1, x0)) - abs(x1 - x0)
          bounded = TRUE
        }
        x[i+1] = x2

      }else if(i == 3){
        posMin = which.min(abs(target-y))
        cat("Best",name.par,"found ",x[posMin],"--",name, "=", y[posMin], "\n", sep=" ")
        
        output = list(
          objective.function = y[posMin],
          par.value = x[posMin]
        )
        return(output)
      }
    }
    
    if(y[i] <= target){
      x[i+1] = x[i] * 1.05
    }else{
      x[i+1] = x[i] * 0.95
    }
  }
  
}

d = calibrateARL(targetARL=168, targetmRL=NULL,n, m, theta=NULL, Ftheta=NULL, dist, mu, sigma, dist.par=dist.par, initial.par = chart.par, replicates=replicates, chart=chart)

