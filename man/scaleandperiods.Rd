\name{scaleandperiods}
\alias{scaleandperiods}
\title{
Scaling data function.
}
\description{
This function scales the time series data in the interval [0,1] as deemed necessary in Nieto-Barajas and Contreras-Cristan (2014) for the time series clustering algorithm. It also obtains the time periods of the data set provided. 
}
\usage{
scaleandperiods(data,scale)
}
\arguments{
  \item{data}{
Data frame with the time series information.
}
  \item{scale}{
Flag that indicates if the time series data should be scaled to the [0,1] interval with a linear transformation as proposed by Nieto-Barajas and Contreras-Cristan (2014). If TRUE, then the time series are scaled to the [0,1] interval. Its value comes directly from the "scale" argument of the clustering functions. 
}
}
\details{
The function considers that the time periods of the data appear as row
names.
}
\value{
\item{periods}{array with the time periods of the data.}
\item{mydata}{data frame with the time series data scaled in [0,1].}
\item{cts}{variable that indicates if some time series were removed because they were constant in time. If no time series were removed, cts = 0. If there were time series removed, cts indicates the column of such time series.} 
}
\author{
Martell-Juarez, D.A.
}
\note{
For internal use.
}