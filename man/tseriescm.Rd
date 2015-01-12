\name{tseriescm}
\alias{tseriescm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for monthly time series clustering. 
}
\description{
Function that perform the time series clustering algorithm described in Nieto-Barajas and Contreras-Cristan (2014) for monthly time series data.
}
\usage{
tseriescm(filedir, maxiter = 2000, p = 1, d = 13, c0eps = 2, c1eps = 1,
          c0beta = 2,c1beta = 1, c0alpha = 2, c1alpha = 1, priora = 0,
          priorb = 0, a = 0.25, b = 0, indlpml = 0)
}
\arguments{
  \item{filedir}{
String with the file directory from which the data will be read. 
}
  \item{maxiter}{
Maximum number of iterations for Gibbs sampling.
}
  \item{p}{
Number of parameters not considered for clustering.
}
  \item{d}{
Number of parameters considered for clustering.
}
  \item{c0eps}{
Shape parameter of the hyper-prior distribution on sig2eps.
}
  \item{c1eps}{
Rate parameter of the hyper-prior distribution on sig2eps. 
}
  \item{c0beta}{
Shape parameter of the hyper-prior distribution on sig2beta. 
}
  \item{c1beta}{
Rate parameter of the hyper-prior distribution on sig2beta. 
}
  \item{c0alpha}{
Shape parameter of the hyper-prior distribution on sig2alpha. 
}
  \item{c1alpha}{
Rate parameter of the hyper-prior distribution on sig2alpha. 
}
  \item{priora}{
Variable that indicates if a prior on parameter "a" is to be assigned. If priora = 0, no prior is assigned. If priora = 1, a prior on "a" is assigned and the user is asked to enter the relevant values for the prior distribution when the function is run.
}
  \item{priorb}{
Variable that indicates if a prior on parameter "b" is to be assigned. If priorb = 0, no prior is assigned.If priorb = 1, a prior on "b" is assigned and the user is asked to enter the relevant values for the prior distribution when the function is run.
}
  \item{a}{
Initial/fixed value of parameter "a".
}
  \item{b}{
Initial/fixed value of parameter "b".
}
  \item{indlpml}{
Variable that indicates if the LPML of the model is to be calculated. If indlpml = 0, LPML is not calculated. If indlpml = 1, LPML is calculated. 
}
}
\details{
It is assumed that the file is saved in csv format. The time periods must be in the first column of the file and there must be a title for each time series in the first row of the file. This function is for internal use. 
}
\value{
\item{mstar}{Number of groups of the chosen cluster configuration.}
\item{gnstar}{Array that contains the group number to which each time series belongs.}
\item{HM}{Heterogeneity Measure of the chosen cluster configuration.} 
\item{arrho}{Acceptance rate of the parameter "rho".}
\item{ara}{Acceptance rate of the parameter "a".}
\item{arb}{Acceptance rate of the parameter "b".}
\item{sig2epssample}{Matrix that in its columns contains the sample of each sig2eps_i's posterior distribution after Gibbs sampling.}
\item{sig2alphasample}{Matrix that in its columns contains the sample of each sig2alpha_i's posterior distribution after Gibbs sampling.}
\item{sig2betasample}{Matrix that in its columns contains the sample of each sig2beta_i's posterior distribution after Gibbs sampling.}
\item{sig2thesample}{Vector that contains the sample of sig2the's posterior distribution after Gibbs sampling.} 
\item{rhosample}{Vector that contains the sample of rho's posterior distribution after Gibbs sampling.}
\item{asample}{Vector that contains the sample of a's posterior distribution after Gibbs sampling.}
\item{bsample}{Vector that contains the sample of b's posterior distribution after Gibbs sampling.}
\item{msample}{Vector that contains the sample of the number of groups at each Gibbs sampling iteration.}
\item{lpml}{If indlpml = 1, lpml contains the value of the LPML of the chosen model.}
}

\examples{
## Do not run
#
# data(stocks)
# write.csv(stocks,file="stocks.csv")
# tseriescm.out <- tseriescm("stocks.csv",maxiter=4000,priora=1,priorb=1)
# 0.5
# 1
# 1
# 1
# 1
#
# The output of the above example is already in the package.
# To show the results:

maxiter <- 4000
data(tseriescm.out)
attach(tseriescm.out)
data(stocks)
write.csv(stocks,file="stocks.csv")
mydata1 <- readscaleperiods("stocks.csv")
mydata <- mydata1$mydata
periods <- mydata1$periods

###### CONVERGENCE ANALYSIS #####

# Trace plots

par(mfrow = c(2,2))
plot(seq((maxiter-1000)/5),tseriescm.out$sig2epssample[,1],type = "l"
,main = "Trace plot of sig2eps",xlab = "Iteration number",ylab = 
"Simulated value")
plot(seq((maxiter-1000)/5),tseriescm.out$sig2alphasample[,1],type = "l"
,main = "Trace plot of sig2alpha",xlab = "Iteration number",ylab = 
"Simulated value")
plot(seq((maxiter-1000)/5),tseriescm.out$sig2betasample[,1],type = "l"
,main = "Trace plot of sig2beta",xlab = "Iteration number",ylab = 
"Simulated value")
plot(seq((maxiter-1000)/5),tseriescm.out$sig2thesample,type = "l",main = 
"Trace plot of sig2theta",xlab = "Iteration number",ylab = "Simulated 
value")

par(mfrow = c(2,2))
plot(seq((maxiter-1000)/5),tseriescm.out$rhosample,type = "l",main = 
"Trace plot of rho",xlab = "Iteration number",ylab = "Simulated value")
plot(seq((maxiter-1000)/5),tseriescm.out$asample,type = "l",main = "Trace 
plot of a",xlab = "Iteration number",ylab = "Simulated value")
plot(seq((maxiter-1000)/5),tseriescm.out$bsample,type = "l",main = "Trace 
plot of b",xlab = "Iteration number",ylab = "Simulated value")
plot(seq((maxiter-1000)/5),tseriescm.out$msample,type = "l",main = "m"
,xlab = "Iteration number",ylab = "Number of groups")

# Histograms

par(mfrow = c(2,2))
hist(tseriescm.out$sig2epssample[,1],main = "Hist. of sig2eps",xlab = 
"Simulated values")
hist(tseriescm.out$sig2alphasample[,1],main = "Hist. of sig2alpha",xlab = 
"Simulated values")
hist(tseriescm.out$sig2betasample[,1],main = "Hist. of sig2beta",xlab = 
"Simulated values")
hist(tseriescm.out$sig2thesample,main = "Hist. of sig2the",xlab = 
"Simulated values")

par(mfrow = c(2,2))
hist(tseriescm.out$rhosample,main = "Hist. of rho",xlab = "Simulated 
values")
hist(tseriescm.out$asample,main = "Hist. of a",xlab = "Simulated values")
hist(tseriescm.out$bsample,main = "Hist. of b",xlab = "Simulated values")
hist(tseriescm.out$msample,main = "Hist. of m",xlab = "Number of groups")

# Ergodic means

par(mfrow = c(2,2))
plot(seq((maxiter-1000)/5),cumsum(tseriescm.out$sig2epssample[,1])/seq
((maxiter-1000)/5),type = "l",main = "Ergodic mean of sig2eps",xlab = 
"Iteration number",ylab = "")
plot(seq((maxiter-1000)/5),cumsum(tseriescm.out$sig2alphasample[,1])/seq
((maxiter-1000)/5),type = "l",main = "Ergodic mean of sig2alpha",xlab = 
"Iteration number",ylab = "")
plot(seq((maxiter-1000)/5),cumsum(tseriescm.out$sig2betasample[,1])/seq
((maxiter-1000)/5),type = "l",main = "Ergodic mean of sig2beta",xlab = 
"Iteration number",ylab = "")
plot(seq((maxiter-1000)/5),cumsum(tseriescm.out$sig2thesample)/seq
((maxiter-1000)/5),type = "l",main = "Ergodic mean of sig2the",xlab = 
"Iteration number",ylab = "")

par(mfrow = c(2,2))
plot(seq((maxiter-1000)/5),cumsum(tseriescm.out$rhosample)/seq((maxiter
-1000)/5),type = "l",main = "Ergodic mean of rho",xlab = "Iteration 
number",ylab = "")
plot(seq((maxiter-1000)/5),cumsum(tseriescm.out$asample)/seq((maxiter
-1000)/5),type = "l",main = "Ergodic mean of a",xlab = "Iteration number"
,ylab = "")
plot(seq((maxiter-1000)/5),cumsum(tseriescm.out$bsample)/seq((maxiter
-1000)/5),type = "l",main = "Ergodic mean of b",xlab = "Iteration number"
,ylab = "")
plot(seq((maxiter-1000)/5),cumsum(tseriescm.out$msample)/seq((maxiter
-1000)/5),type = "l",main = "Ergodic mean of m",xlab = "Iteration number"
,ylab = "")

# Cluster plots 
T <- nrow(mydata)
fT <- floor(T/3)
auxtt <- matrix(0,fT,1)

for(i in 1:fT){
  auxtt[i] <- 3*i
}

for(j in 1:tseriescm.out$mstar){
  
  cc <- as.matrix(which(tseriescm.out$gnstar == j))
  
  cl <- rainbow(nrow(cc))
  
  plot((1:T),mydata[,cc[1,1]],type = "l",main = paste("Group",j),xlab = 
  "",xaxt = 'n',ylab = "Scaled variable in [0,1]",col = cl[1])
  axis(1,at = auxtt,labels = periods[auxtt],las = 2,tck = 0)
  
  if(nrow(cc) > 1){
    for(i in 2:nrow(cc)){
      lines((1:T),mydata[,cc[i,1]],col = cl[i])
    }
  }
  
}

# Console output

cat("Number of groups of the chosen cluster configuration: ",tseriescm.out$mstar,"\n")
for(i in 1:tseriescm.out$mstar){
  cat("Time series in group",i,":",which(tseriescm.out$gnstar == i),"\n")
}
cat("HM Measure: ",tseriescm.out$HM,"\n")
cat("Acceptance rate of rho: ",tseriescm.out$arrho,"\n")
cat("Acceptance rate of a: ",tseriescm.out$ara,"\n")
cat("Acceptance rate of b: ",tseriescm.out$arb,"\n")
}
\references{
Nieto-Barajas, L.E. and Contreras-Cristan, A. (2014) A Bayesian Nonparametric Approach for Time Series Clustering. \emph{Bayesian Analysis} \bold{Vol. 9, No. 1} 147--170.
}
\author{
Martell-Juarez, D.A. and Nieto-Barajas, L.E.
}