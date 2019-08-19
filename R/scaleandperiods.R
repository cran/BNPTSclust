scaleandperiods <-
function(data,scale){
    
  # Function that receives a data frame with the time series data and 
  # scales it in the [0,1] interval, if desired by the user.
  # The function considers that the time periods of the data appear 
  # as row names.
  # 
  # IN:
  #
  # data  <- data frame with the time series information.
  # scale   <- Flag that indicates if the time series data should be scaled to the
  #            [0,1] interval with a linear transformation as proposed by 
  #            Nieto-Barajas and Contreras-Cristan (2014). If TRUE, then the time
  #            series are scaled to the [0,1] interval. Its value comes directly from
  #            the "scale" argument of the clustering functions. 
  #
  # OUT:
  # 
  # periods <- array with the time periods of the data. 
  # mydata  <- data frame with the time series data. 
  # cts     <- variable that indicates if some time series were removed
  #            because they were constant in time. If no time series were
  #            removed, cts = 0. If there were time series removed, cts
  #            indicates the column of such time series. 
    
    n <- nrow(data)
    m <- ncol(data)
    
    periods <- rownames(data)
    
    mydata <- data
    colnames(mydata) <- colnames(data)
    names <- colnames(mydata)
    
    maxima <- matrix(0,m,1)
    minima <- matrix(0,m,1)
    
    for (i in 1:m){
      maxima[i,1] = max(mydata[,i])
      minima[i,1] = min(mydata[,i])  
    }
    
    cts <- which(maxima == minima)
    
    if(length(cts) != 0){
      cat("Removing series ",names[cts]," since they are constant.","\n")
      mydata <- mydata[,-cts]
      n <- nrow(mydata)
      m <- ncol(mydata)
      maxima <- matrix(0,m,1)
      minima <- matrix(0,m,1)
      
      for (i in 1:m){
        maxima[i,1] = max(mydata[,i])
        minima[i,1] = min(mydata[,i])  
      }
      
    }
    
    if(scale==TRUE){
      for (j in 1:m){
        m1 = maxima[j,1] - minima[j,1]
        
        for (k in 1:n){
          mydata[k,j] = 1 + (1/m1)*(mydata[k,j] - maxima[j,1])
        } 
    }
 
      
    }
    
    return(list(periods = periods, mydata = mydata, cts = cts))
    
  }
