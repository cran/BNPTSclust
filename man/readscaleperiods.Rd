\name{readscaleperiods}
\alias{readscaleperiods}
\title{
Scaling data function
}
\description{
This function scales the time series data in the interval [0,1] as deemed necessary in Nieto-Barajas and Contreras-Cristan (2014) for the time series clustering algorithm. It also obtains the time periods of the data set provided. 
}
\usage{
readscaleperiods(filedir)
}
\arguments{
  \item{filedir}{string with the pathfile where the data is saved.}
}
\details{
It is assumed that the file is saved in csv format. The time periods must be in the first column of the file and there must be a title for each time series in the first row of the file. This function is for internal use. 
}
\value{
\item{periods}{array with the time periods of the data.}
\item{mydata}{matrix with the time series data scaled in [0,1].}
\item{cts}{  # cts     <- variable that indicates if some time series were removed because they were constant in time. If no time series were removed, cts = 0. If there were time series removed, cts indicates the column of such time series.}
}
\author{
Martell-Juarez, D.A.
}
\examples{
## The function is currently defined as
readscaleperiods <- function(filedir){
  
  # Function that reads the file with the time series data and scales it
  # in the [0,1] interval.
  # The function considers that the time periods of the data appear on the
  # first column of the file and that the title of each time series 
  # appears on the first row of the file. 
  # The file must be saved in "csv" format. 
  # 
  # IN:
  #
  # filedir <- string with the file directory from which the data will
  #            be read.
  # OUT:
  # 
  # periods <- array with the time periods of the data. 
  # mydata  <- matrix with the time series data. 
  # cts     <- variable that indicates if some time series were removed
  #            because they were constant in time. If no time series were
  #            removed, cts = 0. If there were time series removed, cts
  #            indicates the column of such time series.
  
  mydata <- read.csv(filedir)
  
  periods <- mydata[,1]
  mydata <- data.matrix(mydata[,-1])
  
  n <- nrow(mydata)
  m <- ncol(mydata)
  maxima <- matrix(0,m,1)
  minima <- matrix(0,m,1)
  
  for (i in seq(m)){
    maxima[i,1] = max(mydata[,i])
    minima[i,1] = min(mydata[,i])  
  }
  
  cts <- which(maxima == minima)
  
  if(length(cts) != 0){
    mydata <- mydata[,-cts]
    n <- nrow(mydata)
    m <- ncol(mydata)
    maxima <- matrix(0,m,1)
    minima <- matrix(0,m,1)
    
    for (i in seq(m)){
      maxima[i,1] = max(mydata[,i])
      minima[i,1] = min(mydata[,i])  
    }
    
  }
  
  for (j in seq(m)){
    m = maxima[j,1] - minima[j,1]
    
    for (k in seq(n)){
      mydata[k,j] = 1 + (1/m)*(mydata[k,j] - maxima[j,1])
    }  
      
  }
  
  return(list(periods = periods, mydata = mydata, cts = cts))
  
}
}