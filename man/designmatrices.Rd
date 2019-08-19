\name{designmatrices}
\alias{designmatrices}
\title{
Function that creates the design matrices necessary for the clustering algorithm to work. 
}
\description{
Function that generates the design matrices of the clustering algorithm based on the parameters that the user wants to consider, i.e. level, polinomial trend and/or seasonal components. It also returns the number of parameters that are considered and not considered for clustering. 
}
\usage{
designmatrices(level, trend, seasonality, deg, T, n, fun)
}
\arguments{
  \item{level}{
Variable that indicates if the level of the time series will be considered for clustering. If level = 0, then it is omitted. If level = 1, then it is taken into account.
}
  \item{trend}{
Variable that indicates if the polinomial trend of the model will be considered for clustering. If trend = 0, then it is omitted. If trend = 1, then it is taken into account.
}
  \item{seasonality}{
Variable that indicates if the seasonal components of the model will be considered for clustering. If seasonality = 0, then they are omitted. If seasonality = 1, then they are taken into account.
}
  \item{deg}{
Degree of the polinomial trend of the model.
}
  \item{T}{
Number of periods of the time series.
}
  \item{n}{
Number of time series.
}
  \item{fun}{
Clustering function being used.
}
}
\value{
\item{Z}{Design matrix of the parameters not considered for clustering.}
\item{X}{Design matrix of the parameters considered for clustering.}
\item{p}{Number of parameters not considered for clustering.}
\item{d}{Number of parameters considered for clustering.}
}
\author{
Martell-Juarez, D.A.
}
\note{
For internal use. 
}
