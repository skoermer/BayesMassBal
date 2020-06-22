chib.component <- function(s.l,X,y){

  beta <- do.call(rbind,s.l$beta)
  Sig <- s.l$Sig
  x.unit <- X
  Bprior <- s.l$priors$beta
  Sprior <- s.l$priors$Sig$S
  nu0 <- s.l$priors$Sig$nu0

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


  for(t in 1:T){
    for(i in 1:M){
      s <- Sig[[i]][,t]
      S[upper.tri(S, diag = TRUE)] <- s
      S[lower.tri(S)] <- t(S)[lower.tri(S)]
      Si <- solve(S)
      #Si.full[[i]] <- Si
      xtSix <- t(x.unit) %*% Si %*% x.unit * K
      #xtSixi <- solve(xtSix)

      S0i.use <- diag(S0i)[((i-1)*p.M+1):(i*p.M)]
      S0iB0.use <- S0iB0[((i-1)*p.M+1):(i*p.M)]

      cov.use <- solve(xtSix + diag(S0i.use))
      bhat <- cov.use %*% (S0iB0.use + t(x.unit)%*%Si %*% yibar[[i]] * K)

      p.int[i] <- dtmvnorm(B.bar[((i-1)*p.M + 1):(i*p.M)], mean = bhat[,1], sigma = cov.use, lower = rep(0, times= p.M), log = TRUE)
    }

    postB[t] <- sum(p.int)

  }

  lpostB <- log(mean(exp(postB)))


  for(i in 1:M){
    s <- Sig.bar[[i]]
    S[upper.tri(S, diag = TRUE)] <- s
    S[lower.tri(S)] <- t(S)[lower.tri(S)]

    Sig.barfull[[i]] <- S

    B.use <- B.bar[(((i-1)*p.M)+1):(i*p.M)]
    xB <- x.unit %*% B.use
    xB <- xB[,1]

    psi <- matrix(0, ncol = N, nrow = N)
    for(j in 1:K){
      psi <- (y[[i]][,j] - xB) %*% t(y[[i]][,j] - xB) + psi
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
