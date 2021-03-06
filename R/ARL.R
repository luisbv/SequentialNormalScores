source("NS.r")
source("GetDist.r")

# runs only when script is run by itself
if (sys.nframe() == 0){
  n = 5 #subgroup size
  m = 1000 #reference-sample size
  dist = "Normal"
  mu = c(0, 0) # c(reference sample mean, monitoring sample mean)
  sigma = c(1, 1) # c(reference sample sd, monitoring sample sd)
  
  #### Distribution parameters
  dist.par = c(0, 1, 1) # c(location, scale, shape)
  
  #### Control chart parameters
  #chart.par = c(3); chart = "Shewhart"
  #Shewhart c(k), where k comes from UCL = mu + k*sigma, LCL = mu - k*sigma.
  
  chart.par = c(0.5, 5, 1); chart = "CUSUM"
  #CUSUM c(k, h, t), where k threshold and h is the decision limit. t is the type of the chart (1:positive, 2:negative, 3:two sides)
  
  #chart.par = c(0.2, 2.962); chart = "EWMA"
  #EWMA c(lambda, L), where lambda is the smoothing constant and L multiplies standard deviation to get control limit
  
  #### Other Parameters
  replicates = 20
  print.RL = TRUE
  
  
  a = getARLSNS(n, m, theta=NULL, Ftheta=NULL, dist, mu, sigma, dist.par=dist.par, chart.par=chart.par, print.RL = print.RL, replicates=replicates, chart)
}


getARLSNS <- function(n, m, theta=NULL, Ftheta=NULL, dist, mu, sigma, dist.par=c(0, 1, 1), chart.par, print.RL = FALSE, replicates=10000, chart, progress=TRUE){
  RLs = NULL
  t0 = Sys.time()
  for(r in 1:replicates){
    
    RL = get.RL(n, m, theta, Ftheta, dist, mu, sigma, dist.par, chart=chart, chart.par=chart.par)
    
    RLs = c(RLs, RL)
    
    #print out progress
    if(progress){#if is TRUE
      if(r %% 10 == 0){# every 10 replicates
        t1 = Sys.time()
        remaining.iterations = replicates - r
        remaining.time = remaining.iterations * difftime(t1,t0,units="min") / r
        cat("ARL",round(mean(RLs), digits = 1),"-- SDRL", round(sd(RLs), digits = 1),"--> Time remaining", remaining.time, "in minutes to complete", remaining.iterations, "iterations","\n",sep=" ")
      }
    }
  }
  
  output = list(
    ARL = mean(RLs),
    SDRL = sd(RLs),
    mRL = median(RLs),
    RL.quantiles = quantile(RLs)
  )
  if(print.RL) output$RL = RLs
  
  if(progress) cat("Final ARL", round(mean(RLs), digits = 1), "-- SDRL", round(sd(RLs), digits = 1), "\n", "See output variable for more.\n\n", sep=" ")
  
  return(output)
}


get.RL <- function(n, m, theta=NULL, Ftheta=NULL, dist, mu, sigma, dist.par, chart, chart.par){
  #initilize the reference sample

  Y = NULL
  if(m > 0){#if there are reference sample
    #generate the reference sample
    Y = getDist(n=m, dist = dist, mu=mu[1], stdev=sigma[1], par.location = dist.par[1], par.scale = dist.par[2], par.shape = dist.par[3])
  }
  RL = 0
  in.Control = TRUE
  
  switch(chart,
    Shewhart = {
      k = chart.par[1]
    },
    CUSUM = {
      k = chart.par[1]
      h = chart.par[2]
      type = chart.par[3]
      Cplus = 0
      Cminus = 0
    },
    EWMA = {
      lambda = chart.par[1]
      L = chart.par[2]
      E = 0
    }
  )

  while(in.Control){
    #add one iteration to run length
    RL = RL + 1
    
    #generate the subgroup to monitor
    X = getDist(n=n, dist = dist, mu = mu[2], stdev=sigma[2], par.location = dist.par[1], par.scale = dist.par[2], par.shape = dist.par[3])
    
    #get the normal scores
    Z = NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
    Z = mean(Z)
    
    #if the subgroup is out of the limits
    # an alarm is detected
    switch(chart,
      Shewhart = {
        #if the subgroup is out of the limits an alarm is detected
        if (abs(Z) >= k/sqrt(n)) in.Control = FALSE
      },
      CUSUM = {
        switch(type,
          "1" = {
           Cplus = max(c(0, Cplus + Z*sqrt(n) - k))
          },
          "2" = {
           Cminus = min(c(0, Cminus + Z*sqrt(n) + k))
          },
          "3" = {
           Cplus = max(c(0, Cplus + Z*sqrt(n) - k))
           Cminus = min(c(0, Cminus + Z*sqrt(n) + k))
          }
        )
        
        if (Cplus >= h || Cminus <= -h) in.Control = FALSE
      },
      EWMA = {
        E = lambda*Z + (1-lambda)*E
        
        UCL = L/sqrt(n)*sqrt(lambda/(2-lambda)*(1-(1-lambda)^(2*RL)))
        #LCL = - UCL
        
        if(abs(E) >= UCL) in.Control = FALSE
      }
    )
    
    #update the reference sample
    Y = c(Y, X)
  }
  return(RL)
}

