component_sig <- function(X,y,priors,BTE = c(3000,100000,1), verb = 1){

  ## Tests
  K <- ncol(y[[1]])
  ## Sample Locations
  N <- nrow(y[[1]])
  ## Components
  M <- length(y)

  burn <- BTE[1]
  iters <-BTE[2]
  thin <- BTE[3]

  x.unit <- X

  X.unit <- as.matrix(bdiag(replicate(M,x.unit,simplify = FALSE)))
  X <- do.call("rbind", replicate(K, X.unit, simplify=FALSE))

  no.betas <- ncol(X)/M

  X.M <- list()
  Sig <- list()
  Sig.rows <- length(diag(N)[upper.tri(diag(N),diag = TRUE)])
  for(i in 1:M){
    X.M[[i]] <- X.unit[((i-1)*N+1):(i*N),]
    Sig[[i]] <- matrix(NA, nrow = Sig.rows, ncol = iters)
    Sig[[i]][,1] <- diag(N)[upper.tri(diag(N),diag = TRUE)]
  }

  ybar <- lapply(y,rowMeans)

  Y <- do.call(rbind,y)
  Y <- as.vector(Y)

  beta <- replicate(M, matrix(NA, nrow = no.betas, ncol = iters), simplify = FALSE)

  names(beta) <- names(y)

  bhat <- as.matrix(solve(t(X) %*% X) %*% t(X) %*% Y)
  bhat[which(bhat < 0)] <- 1000

  b.use <- bhat

  if(is.na(priors)){
  S.prior <- lapply(y,function(X){diag((apply(X,1,sd))^2)})
  nu0 <- N
  mu0 <- rep(0, times = no.betas*M)
  V0i <- diag(no.betas*M)/10000000
  }else{
    nu0 <- priors$Sig$nu0
    mu0 <- priors$beta$mu0
    V0i <- solve(priors$beta$V0)
    S.prior <- priors$Sig$S
    }

  V0imu0 <- V0i %*% mu0

  beta.temp <- rep(NA, times = no.betas*M)
  if(verb != 0){
  message("Sampling from posterior distributions")
  pb <- txtProgressBar(min = 0, max = iters/100, initial = 0, style = 3)
  step <- 0
  }

  for(t in 2:iters){
    if(verb != 0 & (t/100) %% 1 == 0){
      step <- step + 1
      setTxtProgressBar(pb,value = step)
    }
    for(i in 1:M){
      xB <- as.vector(X.M[[i]] %*% b.use)
      psi <- matrix(0, ncol = N, nrow = N)
      for(k in 1:K){
        psi <- (y[[i]][,k] - xB) %*% t(y[[i]][,k] - xB) + psi
        }
      psi <- psi + S.prior[[i]]*(nu0)
      W <- rinvwishart(K + nu0, psi)
      Sig[[i]][,t] <- W[upper.tri(W, diag = TRUE)]
      Wi <- solve(W)
      xtWix <- t(x.unit) %*% Wi %*% x.unit * K
      precis <- xtWix + diag(diag(V0i)[((i - 1)*no.betas +1):(i*no.betas)])
      cov.use <- solve(precis)
      bhat <- cov.use %*% (V0imu0[((i - 1)*no.betas +1):(i*no.betas)] + t(x.unit)%*%Wi %*% ybar[[i]] * K)
      beta.temp[((i - 1)*no.betas +1):(i*no.betas)] <- beta[[i]][,t] <- as.vector(
                                         rtmvnorm(n = 1, mean = bhat[,1], sigma = cov.use, lower = rep(0, length = no.betas)))
    }
    b.use <- beta.temp
  }
  if(verb != 0){close(pb)}

  Sig <- lapply(Sig,function(X, b){X[,-(1:b)]}, b = burn)
  Sig <- lapply(Sig, function(X,t){X[,seq(from = 1, to = ncol(X), by = t)]}, t = thin)

  beta <- lapply(beta,function(X, b){X[,-(1:b)]}, b = burn)
  beta <- lapply(beta, function(X,t){X[,seq(from = 1, to = ncol(X), by = t)]}, t = thin)

  return(list(beta = beta,
              Sig = Sig,
              priors= list(beta = list(mu0 = mu0, V0 = solve(V0i)), Sig = list(nu0 = nu0, S = S.prior)), cov.structure = "component"))
}
