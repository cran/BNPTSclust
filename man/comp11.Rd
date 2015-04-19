\name{comp11}
\alias{comp11}
\title{
Univariate ties function
}
\description{
Computes the distinct observations and frequencies in a numeric vector.
}
\usage{
comp11(y)
}
\arguments{
  \item{y}{
  Numeric vector. 
}
}
\details{
The code of the function is the same as the "comp1" function from the
"BNPdensity" package. The change is in the output of the function. This function is for internal use. 
}
\value{
\item{jstar}{variable that rearranges "y" into a vector with its unique values.}
\item{nstar}{frequency of each distinct observation in "y".}
\item{rstar}{number of distinct observations in "y".}
\item{gn}{variable that indicates the group number to which every entry in "y" belongs.}
}
\author{
Martell-Juarez, D.A., Barrios, E., Nieto-Barajas, L. and Pruenster, I.
}
\note{
For internal use.
}