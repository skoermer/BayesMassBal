chib.indep <- function(s.l,X,y, verb){

  beta <- do.call(rbind,s.l$beta)
  Sig <- s.l$Sig
  x.unit <- X
  Bprior <- s.l$priors$beta
  Sprior <- s.l$priors$phi

  ### Components
  M <- length(Sig)

  ### Locations
  N <- nrow(y[[1]])

  ### Number of tests

  K <- ncol(y[[1]])

  T <- ncol(beta)
  P <- nrow(beta)

  p.M <- P/M

  ### Define Ybar and X

  Y <- do.call(rbind,y)
  Ybar <- apply(Y,1,mean)

  Y <- as.vector(Y)

  yibar <- lapply(y,rowMeans)

  X.unit <- bdiag(replicate(M,x.unit,simplify = FALSE))
  X <- do.call("rbind", replicate(K, X.unit, simplify=FALSE))
  ## Prior Hyperparameters

  S0 <- Bprior$V0
  B0 <- Bprior$mu0
  S0i <- solve(S0)
  S0iB0 <- S0i %*% B0

  ### Initialize

  p.int <- rep(NA, times = M)

  B.bar <- apply(beta,1,mean)

  Sig.bar <- list()

  for(j in 1:M){
    Sig.bar[[j]] <- apply(Sig[[j]],1,mean)
  }

  Sig.barfull <- list()

  Si.full <- Bcov.full <- list()
  S <- matrix(NA, N,N)

  postB <- rep(NA, times = T)
  lpost.Sig <- rep(NA, times = M)
  lprior.Sig <- rep(NA, times = M)
  rate <- rep(NA, times = N)

  if(verb != 0){
    message("Approximating integral for log-marignal likelihood")
    pb <- txtProgressBar(min = 0, max = T/100, initial = 0, style = 3)
    step <- 0
  }

  for(t in 1:T){

    if(verb != 0 & (t/100) %% 1 == 0){
      step <- step + 1
      setTxtProgressBar(pb,value = step)
    }
    for(i in 1:M){
      Si <- diag(Sig[[i]][,t])
      xtSix <-  t(x.unit) %*% Si %*% x.unit * K
      #xtSixi <- solve(xtSix)

      S0i.use <- diag(S0i)[((i-1)*p.M+1):(i*p.M)]
      S0iB0.use <- S0iB0[((i-1)*p.M+1):(i*p.M)]

      cov.use <- solve(xtSix + diag(S0i.use))
      bhat <- cov.use %*% (S0iB0.use + t(x.unit)%*%Si %*% yibar[[i]] * K)

      p.int[i] <- dtmvnorm(B.bar[((i-1)*p.M + 1):(i*p.M)], mean = bhat[,1], sigma = cov.use, lower = rep(0, times= p.M), log = TRUE)
    }

    postB[t] <- sum(p.int)

  }
  if(verb != 0){close(pb)}
  lpostB <- log(mean(exp(postB)))


  for(i in 1:M){
    s <- Sig.bar[[i]]

    B.use <- B.bar[(((i-1)*p.M)+1):(i*p.M)]
    xB <- x.unit %*% B.use

    for(j in 1:N){
      rate[j] <- 0.5*(t(xB[j] - y[[i]][j,]) %*% (xB[j] - y[[i]][j,]) )+ Sprior[[i]]$b[j]
    }

    shape <- N/2 + Sprior[[i]]$a

    lpost.Sig[i] <- sum(dgamma(s,shape = shape, rate = rate, log = TRUE))
    lprior.Sig[i] <- sum(dgamma(s,shape = Sprior[[i]]$a,rate = Sprior[[i]]$b, log = TRUE))
  }

  lpriorB <- dtmvnorm(B.bar, mean = B0, sigma = S0, lower =rep(0, times = P), log = TRUE)

  lpost <- sum(lpost.Sig) + lpostB

  lprior <- sum(lprior.Sig) + lpriorB

  Omega <- do.call(c,Sig.bar)
  Omega <- 1/Omega
  Omega <- rep(Omega, times = K)
  Omega <- diag(Omega)

  llik <- dtmvnorm(Y,mean = as.numeric(X%*%B.bar), sigma = Omega, lower =rep(0, times = length(Y)), log = TRUE)

  lml <- llik + lprior - lpost
  return(list(lpost = lpost, llik = llik, lprior =lprior, lml= lml))

}
