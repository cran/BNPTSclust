tseriescq <-
function(filedir,maxiter=2000,p=1,d=5,c0eps=2,c1eps=1,
                      c0beta=2,c1beta=1,c0alpha=2,c1alpha=1,priora=0,
                      priorb=0,a=0.25,b=0,indlpml=0){
  
  # Function that performs the time series clustering algorithm
  # described in Nieto-Barajas and Contreras-Cristan (2014) "A Bayesian
  # Non-Parametric Approach for Time Series Clustering". Bayesian 
  # Analysis, Vol. 9, No. 1 (2014) pp.147-170". This function is 
  # designed for quarterly time series data. 
  #
  # IN:
  #
  # filedir <- string with the file directory from which the data will
  #            be read. 
  # maxiter <- maximum number of iterations for Gibbs sampling.
  #            Default value = 2000.
  # p       <- Number of parameters not considered for clustering.
  #            Default value = 1. 
  # d       <- Number of parameters considered for clustering.
  #            Default value = 5. 
  # c0eps   <- Shape parameter of the hyper-prior distribution 
  #            on sig2eps. Default value = 2.
  # c1eps   <- Rate parameter of the hyper-prior distribution 
  #            on sig2eps. Default value = 1.
  # c0beta  <- Shape parameter of the hyper-prior distribution 
  #            on sig2beta. Default value = 2.
  # c1beta  <- Rate parameter of the hyper-prior distribution 
  #            on sig2beta. Default value = 1.
  # c0alpha <- Shape parameter of the hyper-prior distribution 
  #            on sig2alpha. Default value = 2.
  # c1alpha <- Rate parameter of the hyper-prior distribution 
  #            on sig2alpha. Default value = 1.
  # priora  <- Variable that indicates if a prior on parameter "a" is
  #            to be assigned. If priora = 0, no prior is assigned.
  #            If priora = 1, a prior on "a" is assigned and the user is 
  #            asked to enter the relevant values for the prior 
  #            distribution. Default value = 0.
  # priorb  <- Variable that indicates if a prior on parameter "b" is
  #            to be assigned. If priorb = 0, no prior is assigned.
  #            If priorb = 1, a prior on "b" is assigned and the user is 
  #            asked to enter the relevant values for the prior 
  #            distribution. Default value = 0.
  # a       <- Initial/fixed value of parameter "a". 
  #            Default value = 0.25.
  # b       <- Initial/fixed value of parameter "b". 
  #            Default value = 0.
  # indlpml <- Variable that indicates if the LPML is to be calculated.
  #            If indlpml = 0, LPML is not calculated. If indlpml = 1,
  #            LPML is calculated. Default value = 0.
  # 
  # OUT:
  #
  # mstar   <- Number of groups of the chosen cluster configuration.
  # gnstar  <- Array that contains the group number to which each time
  #            series belongs. 
  # HM      <- Heterogeneity Measure of the chosen cluster configuration.
  # arrho   <- Acceptance rate of the parameter "rho".
  # ara     <- Acceptance rate of the parameter "a".
  # arb     <- Acceptance rate of the parameter "b".
  # sig2epssample <- Matrix that in its columns contains the sample of each sig2eps_i's posterior distribution after Gibbs sampling.
  # sig2alphasample <- Matrix that in its columns contains the sample of each sig2alpha_i's posterior distribution after Gibbs sampling.
  # sig2betasample <- Matrix that in its columns contains the sample of each sig2beta_i's posterior distribution after Gibbs sampling.
  # sig2thesample <- Vector that contains the sample of sig2the's posterior distribution after Gibbs sampling. 
  # rhosample <- Vector that contains the sample of rho's posterior distribution after Gibbs sampling.
  # asample <- Vector that contains the sample of a's posterior distribution after Gibbs sampling.
  # bsample <- Vector that contains the sample of b's posterior distribution after Gibbs sampling.
  # msample <- Vector that contains the sample of the number of groups at each Gibbs sampling iteration.
  # lpml    <- If indlpml = 1, lpml contains the value of the LPML of the
  #            chosen model.
  
  
  if(p%%1 != 0 | d%%1 != 0 | p <= 0 | d < 3){
    stop("p and d must be positive integer numbers and d must be greater than or equal to 3.")
  }
  
  if(maxiter <= 1000){
    stop("maxiter must be greater than 1000.")
  }
  
  if(c0eps <= 0 | c1eps <= 0 | c0beta <= 0 | c1beta <= 0 |
       c0alpha <= 0 | c1alpha <= 0){
    stop("c0eps,c1eps,c0beta,c1beta,c0alpha and c1alpha must be 
         positive numbers.")
  }
  
  if(priora != 0 & priora != 1){
    stop("priora must be either zero (default) or one.")
  }
  
  if(priorb != 0 & priorb != 1){
    stop("priorb must be either zero (default) or one.")
  }
  
  if(indlpml != 0 & indlpml != 1){
    stop("indlpml must be either zero (default) or one.")
  }
  
   
  
  data <- readscaleperiods(filedir)  
  mydata <- data$mydata              # Matrix with the scaled data.
  periods <- data$periods            # Array with the data periods.
  cts <- data$cts                    # Variable that indicates if any time series
                                     # were removed from the original data set because they were constant.
  
  ##### CONSTRUCTION OF THE DESIGN MATRICES #####
  
  T <- nrow(mydata)                        # Number of periods of the time series
  n <- ncol(mydata)                        # Number of time series present in the data
  Z <- matrix(1,T,p)                       # Design matrix for the parameters that will not be used for clustering
  
  if (p > 1){                              # The excluded parameters are the level and the trend terms
    
    for (i in 2:p){
      Z[,i] <- matrix(seq(T)^(i-1),T,1)  
    }
    
    if(d > 3){
      X1 <- matrix(0,T,(d-3))               # Columns for the trend terms in the design matrix. 
      
      for (i in seq((d-3))){
        X1[,i] <- matrix(seq(T)^(i+(p-1)),T,1)  
      }
      
      num <- floor(T/4)                     # Number of years present in the data
      
      if (num < 1){                          # If the number of quarters in the data is less than 4, the design matrix is filled this way
        X2 <- diag(1,(T-1))
        X2 <- cbind(X2,matrix(0,(T-1),1))
        X2 <- rbind(X2,matrix(0,1,T))
        X <- cbind(X1,X2)
      } else{
        X21 <- rbind(diag(1,3),matrix(0,1,3))  # Matrix that contains the indicator functions for the 3 quarters and one row of zeros to avoid singularity problems in the design matrix 
        X2 <- X21
        resid <- T %% 4                         # Number of the quarter (num+1) present in the data
        
        if (num >= 2){
          for (i in 2:num){
            X2 <- rbind(X2,X21)   
          }  
        }
        
        X2 <- rbind(X2,X21[0:resid,])
        X <- cbind(X1,X2)                    # Final presentation of the design matrix of the parameters that will be used for clustering. 
      }
    }else{
      num <- floor(T/4)                     # Number of years present in the data
      
      if (num < 1){                          # If the number of quarters in the data is less than 4, the design matrix is filled this way
        X2 <- diag(1,(T-1))
        X2 <- cbind(X2,matrix(0,(T-1),1))
        X <- cbind(X1,X2)
        X <- rbind(X,matrix(0,1,T))  
        
      } else{
        X21 <- rbind(diag(1,3),matrix(0,1,3))  # Matrix that contains the indicator functions for the 3 quarters and one row of zeros to avoid singularity problems in the design matrix 
        X2 <- X21
        resid <- T %% 4                         # Number of the quarter (num+1) present in the data
        
        if (num >= 2){
          for (i in 2:num){
            X2 <- rbind(X2,X21)   
          }  
        }
        
        X <- rbind(X2,X21[0:resid,])            # Final presentation of the design matrix of the parameters that will be used for clustering. 
      }
    }
  }
  
  if (p == 1){                             # The only excluded parameter is the level of the data
    
    X1 <- matrix(0,T,(d-3))               # Columns for the trend terms in the design matrix. 
    for (i in seq((d-3))){
      X1[,i] <- matrix(seq(T)^(i),T,1)  
    }
    
    num <- floor(T/4)                       # Number of years present in the data
    
    if (num < 1){                            # If the number of quarters in the data is less than 4, the design matrix is filled this way
      X2 <- diag(1,(T-1))
      X2 <- cbind(X2,matrix(0,(T-1),1))
      X2 <- rbind(X2,matrix(0,1,T))  
      X <- cbind(X1,X2)
      
    } else{
      X21 <- rbind(diag(1,3),matrix(0,1,3))  # Matrix that contains the indicator functions for the 3 quarters and one row of zeros to avoid singularity problems in the design matrix 
      X2 <- X21
      resid <- T %% 4                         # Number of the quarter (num+1) present in the data
      
      if (num >= 2){
        for (i in 2:num){
          X2 <- rbind(X2,X21)   
        }  
      }
      
      X2 <- rbind(X2,X21[0:resid,])            # Construction of the quarter indicator matrix according to how many years are present in the data
      
      X <- cbind(X1,X2)                        # Final presentation of the design matrix of the parameters that will be used for clustering
    }
  }
  
  
  ##### SETTING UP VALUES FOR THE PARAMETERS OF THE HYPER-PRIOR DISTRIBUTIONS #####
  
  if(priora == 1){
    print("Enter the value of the mixing proportion of the prior distribution on 'a'")
    pia <- scan("",nmax = 1)
    
    while(pia <= 0 | pia >= 1){
      print("The mixing proportion must be a number in (0,1)")
      pia <- scan("",nmax = 1)
    }
    
    print("Enter the value for the shape parameters of the continuous part of the prior distribution on 'a'")
    q0a <- scan("",nmax = 1)
    q1a <- scan("",nmax = 1)
    
    while(q0a <= 0 | q1a <= 0){
      print("The shape parameters must be positive numbers. Please enter them again.")
      q0a <- scan("",nmax = 1)
      q1a <- scan("",nmax = 1)
    }     
  }
  
  if(priorb == 1){    
    print("Enter the value for the shape parameters of the prior distribution on 'b'")
    q0b <- scan("",nmax = 1)
    q1b <- scan("",nmax = 1)
    
    while(q0b <= 0 | q1b <= 0){
      print("The shape parameters must be positive numbers. Please enter them again.")
      q0b <- scan("",nmax = 1)
      q1b <- scan("",nmax = 1)
    }    
  }
  
  
  ##### INITIAL VALUES FOR THE PARAMETERS THAT WILL BE PART OF THE GIBBS SAMPLING #####
  
  sig2eps <- matrix(1,n,1)                        # Vector that has the diagonal entries of the variance-covariance matrix for every epsilon_i.
  
  sig2beta <- matrix(1,d,1)                       # Vector that has the diagonal entries of the variance-covariance matrix for beta. 
  sigmabeta <- diag(c(sig2beta),d,d)              # Variance-covariance matrix for beta.
  invsigmabeta <- diag(1/c(sig2beta),d,d)         # Inverse variance-covariance matrix for beta.
  sig2alpha <- matrix(1,p,1)                      # Vector that has the diagonal entries of the variance-covariance matrix for alpha.
  sigmaalpha <- diag(c(sig2alpha),p,p)            # Variance-covariance matrix for alpha.
  invsigmaalpha <- diag(1/c(sig2alpha),p,p)       # Inverse variance-covariance matrix for alpha.
  
  sig2the <- 1    # Initial value for sig2the.
  rho <- 0        # Initial value for rho.  
  
  P <- matrix(0,T,T)             # Initial matrix P.                        
  
  for (j in seq(T)){
    for (k in seq(T)){
      P[j,k] <- rho^(abs(j-k))  
    }
  }
  
  R <- sig2the*P                 # Initial matrix R.
  
  
  alpha <- matrix(mvrnorm(n,matrix(0,p,1),sigmaalpha),p,n)  # alpha is a matrix with a vector value of alpha for every time series in its columns. 
  beta <- matrix(mvrnorm(n,matrix(0,d,1),sigmabeta),d,n)    # beta is a matrix with a vector value of beta for every time series in its columns.
  theta <- matrix(mvrnorm(n,matrix(0,T,1),R),T,n)           # theta is a matrix with a vector value of theta for every time series in its columns.
  gamma <- rbind(beta,theta)                                # gamma is the union by rows of the beta and theta matrices.  
  
  
  iter <- 0                                    # Counter for each Gibbs sampling iteration.
  iter1 <- 0                                   # Counter for the number of iterations saved during the Gibbs sampling.
  arrho <- 0                                   # Variable that will contain the acceptance rate of rho in the Metropolis-Hastings step.
  ara <- 0                                     # Variable that will contain the acceptance rate of a in the Metropolis-Hastings step.
  arb <- 0                                     # Variable that will contain the acceptance rate of b in the Metropolis-Hastings step.
  sim <- matrix(0,n,n)                         # Initialization of the similarities matrix.
  memory <- matrix(0,((maxiter-1000)/5)*n,n)   # Matrix that will contain the cluster configuration of every iteration that is saved during the Gibbs sampling.
  memorygn <- matrix(0,((maxiter-1000)/5),n)   # Matrix that will save the group number to which each time series belongs in every iteration saved.  
  sig2epssample <- matrix(0,((maxiter-1000)/5),n)   # Matrix that in its columns will contain the sample of each sig2eps_i's posterior distribution after Gibbs sampling.
  sig2alphasample <- matrix(0,((maxiter-1000)/5),p) # Matrix that in its columns will contain the sample of each sig2alpha_i's posterior distribution after Gibbs sampling.
  sig2betasample <- matrix(0,((maxiter-1000)/5),d)  # Matrix that in its columns will contain the sample of each sig2beta_i's posterior distribution after Gibbs sampling.
  sig2thesample <- matrix(0,((maxiter-1000)/5),1)   # Vector that will contain the sample of sig2the's posterior distribution after Gibbs sampling. 
  rhosample <- matrix(0,((maxiter-1000)/5),1)       # Vector that will contain the sample of rho's posterior distribution after Gibbs sampling.
  asample <- matrix(0,((maxiter-1000)/5),1)         # Vector that will contain the sample of a's posterior distribution after Gibbs sampling.
  bsample <- matrix(0,((maxiter-1000)/5),1)         # Vector that will contain the sample of b's posterior distribution after Gibbs sampling.
  msample <- matrix(0,((maxiter-1000)/5),1)         # Vector that will contain the sample of the number of groups at each Gibbs sampling iteration.
  
  if(indlpml != 0){
    iter2 <- 0
    auxlpml <- matrix(0,((maxiter-1000)/10),n)
  }
  
  ##### BEGINNING OF GIBBS SAMPLING #####
  
  while(iter < maxiter){
    
    ##### 1) SIMULATION OF ALPHA'S POSTERIOR DISTRIBUTION #####
    
    for(i in 1:n){
      
      sigmaeps <- diag(c(sig2eps[i]),T)
      Q <- sigmaeps + R
      Qinv <- chol2inv(chol(Q))
      Vinv <- t(X) %*% Qinv %*% X + invsigmabeta
      V <- chol2inv(chol(Vinv))
      
      Winv <- Qinv + (Qinv %*% X %*% V %*% t(X) %*% Qinv)
      W <- chol2inv(chol(Winv))
      
      Valphainv <- (t(Z) %*% Winv %*% Z) + invsigmaalpha
      Valpha <- chol2inv(chol(Valphainv))
      
      mualpha <- Valpha %*% t(Z) %*% Winv %*% mydata[,i]
      
      alpha[,i] <- mvrnorm(1,mualpha,Valpha)
      
    }
    
    
    ##### 2) SIMULATION OF GAMMA'S = (BETA,THETA) POSTERIOR DISTRIBUTION #####
    
    
    for(i in 1:n){
      
      gr <- comp11(gamma[1,-i])                         # Only the first entries of gamma[,-i] are compared to determine the cluster configuration
      jstar <- gr$jstar                                 # Object that contains the positions of the unique vectors in gamma[,-i]
      gmi <- gamma[,-i]                                 # Matrix with all the elements of gamma, except for the i-th element
      gammastar <- as.matrix(gmi[,jstar])               # Matrix with the unique vectors in gamma(-i)
      mi <- gr$rstar                                    # Number of unique vectors in gamma(-i) (Number of groups)
      nstar <- gr$nstar                                 # Frequency of each unique vector in gamma(-i)
      betastar <- as.matrix(gammastar[1:d,])            # Separation of unique vectors between betastar and thetastar
      thetastar <- as.matrix(gammastar[(d+1):(T+d),])
      
      # Matrices necessary for the following steps
      sigmaeps <- sig2eps[i]*diag(1,T)
      invsigmaeps <- (1/sig2eps[i])*diag(1,T)
      
      Q <- sigmaeps + R
      Qinv <- chol2inv(chol(Q))
      Vinv <- t(X) %*% Qinv %*% X + invsigmabeta
      V <- chol2inv(chol(Vinv))
      
      Winv <- Qinv + (Qinv %*% X %*% V %*% t(X) %*% Qinv)
      W <- chol2inv(chol(Winv))
      
      # Computing weigths for gamma(i)'s posterior distribution
      dj <- matrix(0,mi,1)
      d0 <- (b + a*mi)*dmvnorm(mydata[,i],(Z %*% alpha[,i]),W)
      
      den <- 0
      
      for(j in 1:mi){
        dj[j] <- (nstar[j] - a)*dmvnorm(mydata[,i],(Z %*% alpha[,i] + X %*% betastar[,j] + thetastar[,j]),sigmaeps)
      }
      
      den <- d0 + sum(dj)
      if(den == 0){
        d0 <- (b + a*mi)*dmvnorm(mydata[,i],(Z %*% alpha[,i]),W,log=TRUE)
        for(j in 1:mi){
          dj[j] <- (nstar[j] - a)*dmvnorm(mydata[,i],(Z %*% alpha[,i] + X %*% betastar[,j] + thetastar[,j]),sigmaeps,log=TRUE)          
        }
        dj <- rbind(dj,d0)
        aa <- min(dj)
        q <- exp(dj-aa)/sum(exp(dj-aa))
      }else{
        q <- dj/den
        q <- rbind(q,d0/den)
      }
      
      # Sampling a number between 1 and (mi+1) to determine what will be the simulated value for gamma(i)
      # The probabilities of the sample are based on the weights previously computed
      
      y <- sample((1:(mi+1)), size=1, prob = q)
      
      # If sample returns the value (mi+1), a new vector from g0 will be simulated and assigned to gamma(i)
      
      if (y == mi+1){
        Sthetai <- chol2inv(chol(invsigmaeps + chol2inv(chol(R))))
        muthetai <- Sthetai %*% invsigmaeps %*% (mydata[,i] - (Z %*% alpha[,i]) - (X %*% beta[,i]))
        mubetai <- V %*% t(X) %*% Qinv %*% (mydata[,i] - (Z %*% alpha[,i]))
        beta0 <- matrix(mvrnorm(1,mubetai,V),d,1)
        theta0 <- matrix(mvrnorm(1,muthetai,Sthetai),T,1)
        gamma[,i] <- rbind(beta0,theta0)
        
      } else{
        gamma[,i] = gammastar[,y]                     # Otherwise, column y from gammastar will be assigned to gamma(i)
      }
      
    }
    
    
    ##### 2.1) ACCELERATION STEP AND CONSTRUCTION OF SIMILARITIES MATRIX #####
    
    gr <- comp11(gamma[1,])                   # Computation of all latent classes of the gamma vectors after the simulation of their posterior distribution.
    jstar <- gr$jstar
    gammastar <- as.matrix(gamma[,jstar])     # Unique values of the gamma vectors. 
    m <- gr$rstar                             # Total number of latent classes (groups).
    nstar <- gr$nstar                         # Frequency of each latent class (group).
    gn <- gr$gn                               # Identifier of the group to which each time series belongs.
    beta <- as.matrix(gamma[(1:d),])          # Splitting the gamma vectors in beta and theta. 
    theta <- as.matrix(gamma[((d+1):(T+d)),])
    betastar <- as.matrix(gammastar[(1:d),])
    thetastar <- as.matrix(gammastar[((d+1):(T+d)),])
    
    for(j in 1:m){
      
      cc <- which(gn == j)        # Identifying the cluster configuration of each group.
      aux <- matrix(0,T,T)        # Calculating the necessary matrices for the simulation of the distributions for the acceleration step.
      aux1 <- matrix(0,T,1)
      aux2 <- matrix(0,T,1)
      
      for(i in 1:nstar[j]){
        aux <- aux + diag((1/sig2eps[cc[i]]),T)
        aux1 <- aux1 + (diag((1/sig2eps[cc[i]]),T) %*% (mydata[,i] - Z %*% alpha[,i] - X %*% betastar[,j]))
        aux2 <- aux2 + (diag((1/sig2eps[cc[i]]),T) %*% (mydata[,i] - Z %*% alpha[,i] - thetastar[,j]))    
      }
      
      Sthetastar <- chol2inv(chol(aux + chol2inv(chol(R))))
      muthetastar <- Sthetastar %*% aux1
      Sbetastar <- chol2inv(chol((t(X) %*% aux %*% X) + invsigmabeta))
      mubetastar <- Sbetastar %*% t(X) %*% aux2
      
      beta[,cc] <- mvrnorm(1,mubetastar,Sbetastar)
      theta[,cc] <- mvrnorm(1,muthetastar,Sthetastar)
      
      # Computation of similarities matrix and saving the cluster configuration of the current iteration.
      if((iter %% 5) == 0 & iter >= 1000){
        
        for(i1 in 1:nstar[j]){
          for(i2 in i1:nstar[j]){
            sim[cc[i1],cc[i2]] <- sim[cc[i1],cc[i2]] + 1
            sim[cc[i2],cc[i1]] <- sim[cc[i2],cc[i1]] + 1
            memory[(cc[i1] + (n*iter1)),cc[i2]] <- memory[(cc[i1] + (n*iter1)),cc[i2]] + 1
            memory[(cc[i2] + (n*iter1)),cc[i1]] <- memory[(cc[i2] + (n*iter1)),cc[i1]] + 1
          }
        }
        
      }
      
    }
    
    gamma <- rbind(beta,theta)                # Obtaining all gamma vectors after the acceleration step.
    gr <- comp11(gamma[1,])                   # Obtaining all the latent classes in the gamma vectors.
    jstar <- gr$jstar
    gammastar <- as.matrix(gamma[,jstar])     # Unique values of the gamma vectors. 
    m <- gr$rstar                             # Number of groups after acceleration step.
    nstar <- gr$nstar                         # Frequency of each group.
    gn <- gr$gn                               # Identifier of the group to which each latent class belongs.
    beta <- as.matrix(gamma[(1:d),])          # Splitting the gamma vectors between beta and theta.
    theta <- as.matrix(gamma[((d+1):(T+d)),])
    betastar <- as.matrix(gammastar[(1:d),])
    thetastar <- as.matrix(gammastar[((d+1):(T+d)),])
    
    
    ##### 3) SIMULATION OF SIG2EPS' POSTERIOR DISTRIBUTION #####
    
    M <- t(mydata - Z%*%alpha - X%*%beta - theta) %*% (mydata - Z%*%alpha - X%*%beta - theta)
    
    sig2eps <- 1/rgamma(n,(c0eps + T/2),(c1eps + diag(M)/2))
    
    
    ##### 4) SIMULATION OF SIMGAALPHA'S POSTERIOR DISTRIBUTION #####
    
    sig2alpha <- 1/rgamma(p,(c0alpha + n/2),(c1alpha + rowSums(alpha^2)))
    
    sigmaalpha <- diag(c(sig2alpha),p,p)
    invsigmaalpha <- diag(1/c(sig2alpha),p,p)
    
    
    ##### 5) SIMULATION OF SIGMABETA'S POSTERIOR DISTRIBUTION #####
    
    sig2beta <- 1/rgamma(d,(c0beta + m/2),(c1beta + colSums(betastar^2)/2))
    
    sigmabeta <- diag(c(sig2beta),d,d)
    invsigmabeta <- diag(1/c(sig2beta),d,d)
    
    
    ##### 6) SIMULATION OF SIG2THE'S POSTERIOR DISTRIBUTION #####
    
    cholP <- chol(P)              # Calculation of the Cholesky factorization of P.
    Pinv <- chol2inv(cholP)       # Obtaining the inverse of P. 
    s1 <- 0
    
    # Calculating the sum necessary for the rate parameter of the posterior distribution.
    for(j in 1:m){
      s1 <- s1 + t(as.matrix(thetastar[,j])) %*% Pinv %*% as.matrix(thetastar[,j])
    }
    
    sig2the <- 1/rgamma(1,(m*T/2),(s1/2))
    
    
    ##### 7) SIMULATION OF RHO'S POSTERIOR DISTRIBUTION (Metropolis-Hastings step) #####
    
    rhomh <- runif(1,-1,1)            # Sampling from the proposal distribution.
    
    Pmh <- matrix(0,T,T)
    
    # Calculating the matrix P for the proposed value rhomh.
    for (j in 1:T){
      for (k in 1:T){
        Pmh[j,k] <- rhomh^(abs(j-k))  
      }
    }
    
    cholPmh <- chol(Pmh)              # Calculating the Cholesky factor of Pmh.
    Pmhinv <- chol2inv(cholPmh)       # Obtaining the inverse from Pmh
    s <- 0
    
    # Calculating the sum necessary for the computation of the acceptance probability.
    for(j in 1:m){
      s <- s + t(as.matrix(thetastar[,j])) %*% (Pmhinv-Pmh) %*% as.matrix(thetastar[,j])
    }
    
    # Computation of the acceptance probability.
    q <- (-m)*(log(prod(diag(cholPmh)))- log(prod(diag(cholP)))) - ((1/(2*sig2the))*s) + (1/2)*(log(1 + rhomh*rhomh) - log(1 + rho*rho)) - log(1 - rhomh*rhomh) + log(1 - rho*rho) 
    
    # Definition of the acceptance probability. 
    quot <- min(0,q)
    
    # Sampling a uniform random variable in [0,1] to determine if the proposal is accepted or not.
    unif1 <- runif(1,0,1)
    
    # Acceptance step.
    if(log(unif1) <= quot){
      rho <- rhomh
      arrho <- arrho + 1
      
      for (j in seq(T)){
        for (k in seq(T)){
          P[j,k] <- rho^(abs(j-k))  
        }
      }
      
    }
    
    R <- sig2the*P  
    
    
    ##### 8) SIMULATION OF A'S POSTERIOR DISTRIBUTION (METROPOLIS-HASTINGS WITH UNIFORM PROPOSALS) #####
    
    if(priora == 1){
      
      if (b < 0){
        amh <- runif(1,-b,1)  
      } else{
        unif2 <- runif(1,0,1)
        if (unif2 <= 0.5){
          amh <- 0
        } else{
          amh <- runif(1,0,1) 
        }
      }
      
      # If b is not greater than -a, then accept the proposal directly.
      if ((a+b) <= 0){
        a <- amh
        print("a+b < 0")
      } else{
        
        quot1 <- 0
        
        if(m > 1){
          for (j in 1:(m-1)){
            quot1 <- quot1 + log(b + j*amh) + log(gamma(nstar[j] - amh)) - log(gamma(1 - amh)) - log(b + j*a) - log(gamma(nstar[j] - a)) + log(gamma(1 - a))
          }
        }
        
        quot1 <- quot1 + log(gamma((nstar[m] - amh))) - log(gamma(1 - amh)) - log(gamma((nstar[m] - a))) + log(gamma(1 - a))
        
        if (a == 0){
          fa <- 0.5 
        } else{
          fa <- 0.5*dbeta(a,q0a,q1a) 
        }
        
        if (amh == 0){
          famh <- 0.5
        } else{
          famh <- 0.5*dbeta(amh,q0a,q1a)
        }
        
        # Quotient to evaluate the Metropolis-Hastings step in logs    
        quot1 <- quot1 + log(famh) - log(fa)
        
        # Determination of the probability for the Metropolis-Hastings step
        alphamh1 <- min(quot1,0)
        
        unif3 <- runif(1,0,1)
        
        # Acceptance step
        if (log(unif3) <= alphamh1){
          a <- amh
          ara <- ara + 1
        }
        
      }
      
    }
    
    
    ##### 9) SIMULATION OF B'S POSTERIOR DISTRIBUTION (METROPOLIS-HASTINGS WITH GAMMA PROPOSALS) #####
    
    if(priorb == 1){
      
      y1 <- rgamma(1,1,0.1)
      bmh <- y1 - a
      
      # If b is not greater than -a, then accept the proposal directly.
      if ((a+b) <= 0){
        b <- bmh  
        print("a+b < 0")
      } else{
        
        quot2 <- 0
        
        if(m > 1){
          for (j in 1:(m-1)){
            quot2 <- quot2 + log(bmh + j*a) - log(b + j*a)
          }
        }
        
        fb <- dgamma(a+b,q0b,q1b)      
        fbmh <- dgamma(y1,q0b,q1b)
        
        # Quotient to evaluate the Metropolis-Hastings step in logs
        quot2 <- quot2 + (log(gamma(bmh+1)) - log(gamma(bmh+n)) - log(gamma(b+1)) + log(gamma(b+n))) + (log(fbmh) - log(fb)) - 0.1*(b - bmh)
        
        # Determination of the probability for the Metropolis-Hastings step
        alphamh2 <- min(quot2,0)
        
        unif4 <- runif(1,0,1)
        
        # Acceptance step
        if (log(unif4) <= alphamh2){
          b <- bmh
          arb <- arb + 1
        }
        
      }
      
    }
    
    
    # 1000 iterations are left as burn-in period and then 1 in every 5
    # iterations is saved as a sample of the posterior distribution of
    # each parameter
    
    if((iter %% 5) == 0 & iter >= 1000){
      iter1 <- iter1 + 1
      sig2epssample[iter1,] <- sig2eps
      sig2alphasample[iter1,] <- sig2alpha
      sig2betasample[iter1,] <- sig2beta
      sig2thesample[iter1] <- sig2the
      rhosample[iter1] <- rho
      asample[iter1] <- a
      bsample[iter1] <- b
      msample[iter1] <- m
      memorygn[iter1,] <- gn
    }
    
    if(indlpml != 0){
      if((iter %% 10) == 0 & iter >= 1000){
        iter2 <- iter2 + 1
        for(i in 1:n){
          for(j in 1:m){        
            auxlpml[iter2,i] <- auxlpml[iter2,i] + ((nstar[j]-a)/(b+n))*dmvnorm(mydata[,i],(Z %*% alpha[,i] + X %*% betastar[,j] + thetastar[,j]),diag(sig2eps[i],T))
          }
          
          sigmaeps <- diag(sig2eps[i],T)
          invsigmaeps <- diag((1/sig2eps[i]),T)
          
          Q <- sigmaeps + R
          Qinv <- solve(Q)
          Vinv <- t(X) %*% Qinv %*% X + invsigmabeta
          V <- solve(Vinv)
          
          Winv <- Qinv + (Qinv %*% X %*% V %*% t(X) %*% Qinv)
          W <- solve(Winv)
          
          auxlpml[iter2,i] <- auxlpml[iter2,i] + ((b+(a*m))/(b+n))*dmvnorm(mydata[,i],(Z %*% alpha[,i]),W) 
        }
      }
    }
    
    iter <- iter + 1
    if(iter %% 50 == 0){
      cat("Iteration Number: ",iter,". Progress: ",(iter/maxiter)*100,"%","\n") 
    }
    
  }
  
  ##### END OF GIBBS SAMPLING #####
  
  # Calculation of acceptance rates and similarities matrix
  arrho <- arrho/iter
  ara <- ara/iter
  arb <- arb/iter
  sim <- sim/iter1
  
  
  dist <- matrix(0,((maxiter-1000)/5),1)
  
  # Calculating the distance between each cluster configuration to the similarities matrix
  for (i in seq(((maxiter-1000)/5))){
    aux4 <- memory[(((i-1)*n)+1):(i*n),] - sim
    dist[i] <- norm(aux4,"F")
  }
  
  # Determining which cluster configuration minimizes the distance to the similarities matrix
  mstar <- msample[which.min(dist)]
  gnstar <- memorygn[which.min(dist),]
  
  
  ###### CONVERGENCE ANALYSIS #####
  
  # Trace plots
  
  par(mfrow = c(2,2))
  plot(seq((maxiter-1000)/5),sig2epssample[,1],type = "l",main = "Trace plot of sig2eps",xlab = "Iteration number",ylab = "Simulated value")
  plot(seq((maxiter-1000)/5),sig2alphasample[,1],type = "l",main = "Trace plot of sig2alpha",xlab = "Iteration number",ylab = "Simulated value")
  plot(seq((maxiter-1000)/5),sig2betasample[,1],type = "l",main = "Trace plot of sig2beta",xlab = "Iteration number",ylab = "Simulated value")
  plot(seq((maxiter-1000)/5),sig2thesample,type = "l",main = "Trace plot of sig2theta",xlab = "Iteration number",ylab = "Simulated value")
  
  par(mfrow = c(2,2))
  plot(seq((maxiter-1000)/5),rhosample,type = "l",main = "Trace plot of rho",xlab = "Iteration number",ylab = "Simulated value")
  plot(seq((maxiter-1000)/5),asample,type = "l",main = "Trace plot of a",xlab = "Iteration number",ylab = "Simulated value")
  plot(seq((maxiter-1000)/5),bsample,type = "l",main = "Trace plot of b",xlab = "Iteration number",ylab = "Simulated value")
  plot(seq((maxiter-1000)/5),msample,type = "l",main = "m",xlab = "Iteration number",ylab = "Number of groups")
  
  # Histograms
  
  par(mfrow = c(2,2))
  hist(sig2epssample[,1],main = "Hist. of sig2eps",xlab = "Simulated values")
  hist(sig2alphasample[,1],main = "Hist. of sig2alpha",xlab = "Simulated values")
  hist(sig2betasample[,1],main = "Hist. of sig2beta",xlab = "Simulated values")
  hist(sig2thesample,main = "Hist. of sig2the",xlab = "Simulated values")
  
  par(mfrow = c(2,2))
  hist(rhosample,main = "Hist. of rho",xlab = "Simulated values")
  hist(asample,main = "Hist. of a",xlab = "Simulated values")
  hist(bsample,main = "Hist. of b",xlab = "Simulated values")
  hist(msample,main = "Hist. of m",xlab = "Number of groups")
  
  # Ergodic means
  
  par(mfrow = c(2,2))
  plot(seq((maxiter-1000)/5),cumsum(sig2epssample[,1])/seq((maxiter-1000)/5),type = "l",main = "Ergodic mean of sig2eps",xlab = "Iteration number",ylab = "")
  plot(seq((maxiter-1000)/5),cumsum(sig2alphasample[,1])/seq((maxiter-1000)/5),type = "l",main = "Ergodic mean of sig2alpha",xlab = "Iteration number",ylab = "")
  plot(seq((maxiter-1000)/5),cumsum(sig2betasample[,1])/seq((maxiter-1000)/5),type = "l",main = "Ergodic mean of sig2beta",xlab = "Iteration number",ylab = "")
  plot(seq((maxiter-1000)/5),cumsum(sig2thesample)/seq((maxiter-1000)/5),type = "l",main = "Ergodic mean of sig2the",xlab = "Iteration number",ylab = "")
  
  par(mfrow = c(2,2))
  plot(seq((maxiter-1000)/5),cumsum(rhosample)/seq((maxiter-1000)/5),type = "l",main = "Ergodic mean of rho",xlab = "Iteration number",ylab = "")
  plot(seq((maxiter-1000)/5),cumsum(asample)/seq((maxiter-1000)/5),type = "l",main = "Ergodic mean of a",xlab = "Iteration number",ylab = "")
  plot(seq((maxiter-1000)/5),cumsum(bsample)/seq((maxiter-1000)/5),type = "l",main = "Ergodic mean of b",xlab = "Iteration number",ylab = "")
  plot(seq((maxiter-1000)/5),cumsum(msample)/seq((maxiter-1000)/5),type = "l",main = "Ergodic mean of m",xlab = "Iteration number",ylab = "")
  
  
  ##### PLOTTING CLUSTERS AND HM MEASURE CALCULATION #####
  
  HM <- 0
  fT <- floor(T/3)
  auxtt <- matrix(0,fT,1)
  
  for(i in 1:fT){
    auxtt[i] <- 3*i
  }
  
  for(j in 1:mstar){
    
    cc <- as.matrix(which(gnstar == j))
    
    cl <- rainbow(nrow(cc))
    
    plot((1:T),mydata[,cc[1,1]],type = "l",main = paste("Group",j),xlab = "",xaxt = 'n',ylab = "Scaled variable in [0,1]",col = cl[1])
    axis(1,at = auxtt,labels = periods[auxtt],las = 2,tck = 0)
    
    if(nrow(cc) > 1){
      
      for(i in 2:nrow(cc)){
        lines((1:T),mydata[,cc[i,1]],col = cl[i])
      }
      
    }
    
    HM1 <- 0
    
    if(length(cc) > 1){
      for(i1 in 1:length(cc)){
        for(i2 in 1:i1){
          HM1 <- HM1 + sum((mydata[,cc[i2]] - mydata[,cc[i1]])^2)
        }
      }
      
      HM <- HM + (2/(length(cc)-1))*HM1
    }
    
  }
  
  cat("Number of groups of the chosen cluster configuration: ",mstar,"\n")
  for(i in 1:mstar){
    cat("Time series in group",i,":",which(gnstar == i),"\n")
  }
  cat("HM Measure: ",HM,"\n")
  cat("Acceptance rate of rho: ",arrho,"\n")
  cat("Acceptance rate of a: ",ara,"\n")
  cat("Acceptance rate of b: ",arb,"\n")
  if(length(cts) != 0){
    cat("The following time series were removed because they were constant through time: ", cts)
  }
  
  if(indlpml != 0){    
    auxlpml <- 1/auxlpml
    cpo <- colMeans(auxlpml)
    cpo <- 1/cpo
    lpml <- sum(log(cpo))
    cat("LPML: ",lpml,"\n")
  }
  
  if(indlpml !=0){
    return(list(mstar = mstar,gnstar = gnstar,HM = HM,arrho = arrho,ara = ara,arb = arb,lpml = lpml,sig2epssample = sig2epssample,sig2alphasample = sig2alphasample,sig2betasample = sig2betasample,sig2thesample = sig2thesample,rhosample = rhosample,asample = asample,bsample = bsample,msample = msample))
  }else{
    return(list(mstar = mstar,gnstar = gnstar,HM = HM,arrho = arrho,ara = ara,arb = arb,sig2epssample = sig2epssample,sig2alphasample = sig2alphasample,sig2betasample = sig2betasample,sig2thesample = sig2thesample,rhosample = rhosample,asample = asample,bsample = bsample,msample = msample))
  }
  
}
