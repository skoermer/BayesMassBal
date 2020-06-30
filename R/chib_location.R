chib.location <- function(s.l,X,y, verb){

  beta <- do.call(rbind,s.l$beta)
  Sig <- s.l$Sig
  x.unit <- X
  Bprior <- s.l$priors$beta
  Sprior <- s.l$priors$Sig$S
  nu0 <- s.l$priors$Sig$nu0

  ### Components
  M <- length(y)

  ### Locations
  N <- nrow(y[[1]])

  ### Number of tests

  K <- ncol(y[[1]])

  T <- ncol(beta)
  P <- nrow(beta)

  p.M <- P/M



  X.unit <- bdiag(replicate(M,x.unit,simplify = FALSE))
  X.unit <- as.matrix(X.unit)

  X.shuffle <- matrix(1:(N*M), ncol = N, nrow = M, byrow = TRUE)
  X.shuffle <- as.vector(X.shuffle)

  X.unit <- X.unit[X.shuffle,]

  X <- do.call("rbind", replicate(K, X.unit, simplify=FALSE))

  y.reorg <- list()
  X.N <- list()

  for(i in 1:N){
    y.reorg[[i]] <- matrix(NA,ncol = K, nrow = M)
    X.N[[i]] <- X.unit[((i-1)*M+1):(i*M),]
    for(j in 1:M){
      y.reorg[[i]][j,] <- y[[j]][i,]
    }
  }

  ybar <- lapply(y.reorg,rowMeans)
  ybar <- do.call(rbind, ybar)
  ybar <- as.vector(t(ybar))

  Y <- do.call(rbind,y.reorg)
  Y <- as.vector(Y)

  ## Prior Hyperparameters

  S0 <- Bprior$V0
  B0 <- Bprior$mu0
  S0i <- solve(S0)
  S0iB0 <- S0i %*% B0

  ### Initialize



  B.bar <- apply(beta,1,mean)

  Sig.bar <- list()

  for(j in 1:N){
    Sig.bar[[j]] <- apply(Sig[[j]],1,mean)
  }

  Sig.barfull <- list()

  Si.full <- Bcov.full <- list()
  S <- matrix(NA, M,M)

  postB <- rep(NA, times = T)
  lpost.Sig <- rep(NA, times = M)
  lprior.Sig <- rep(NA, times = M)
  Si <- list()

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

    for(i in 1:N){
      s <- Sig[[i]][,t]
      S[upper.tri(S, diag = TRUE)] <- s
      S[lower.tri(S)] <- t(S)[lower.tri(S)]
      Si[[i]] <- solve(S)
    }

    Wi <- as.matrix(bdiag(Si))
    XtWiX <- t(X.unit) %*% Wi %*% X.unit * K

    cov.use <- solve(XtWiX + S0i)
    bhat <- cov.use %*% (S0iB0 + t(X.unit)%*%Wi %*% ybar * K)

    postB[t] <- dtmvnorm(B.bar, mean = bhat[,1], sigma = cov.use, lower = rep(0, times= P))

  }
  if(verb != 0){close(pb)}
  lpostB <- log(mean(postB))


  for(i in 1:N){
    s <- Sig.bar[[i]]
    S[upper.tri(S, diag = TRUE)] <- s
    S[lower.tri(S)] <- t(S)[lower.tri(S)]

    Sig.barfull[[i]] <- S

    B.use <- B.bar[(((i-1)*p.M)+1):(i*p.M)]
    xB <- X.N[[i]] %*% B.bar
    xB <- xB[,1]

    psi <- matrix(0, ncol = M, nrow = M)
    for(j in 1:K){
      psi <- (y.reorg[[i]][,j] - xB) %*% t(y.reorg[[i]][,j] - xB) + psi
    }

    psi <- psi + Sprior[[i]]*(nu0)

    lpost.Sig[i] <- dinvwishart(S,nu0+K,psi, log= TRUE)
    lprior.Sig[i] <- dinvwishart(S,nu0,Sprior[[i]]*nu0, log = TRUE)
  }

  lpriorB <- dtmvnorm(B.bar, mean = B0, sigma = S0, lower =rep(0, times = P), log = TRUE)


  lpost <- sum(lpost.Sig) + lpostB

  lprior <- sum(lprior.Sig) + lpriorB

  Omega <- bdiag(rep(Sig.barfull, times = K))

  llik <- dtmvnorm(Y,mean = as.numeric(X%*%B.bar), sigma = Omega, lower =rep(0, times = length(Y)), log = TRUE)

  lml <- llik + lprior - lpost
  return(list(lpost = lpost, llik = llik, lprior =lprior, lml= lml))

}
