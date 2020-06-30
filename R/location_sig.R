location.sig <- function(X,y,priors,BTE = c(3000,100000,1)){

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

  X.unit <- bdiag(replicate(M,x.unit,simplify = FALSE))
  X.unit <- as.matrix(X.unit)

  X.shuffle <- matrix(1:(N*M), ncol = N, nrow = M, byrow = TRUE)
  X.shuffle <- as.vector(X.shuffle)

  X.unit <- X.unit[X.shuffle,]

  X <- do.call("rbind", replicate(K, X.unit, simplify=FALSE))

  no.betas <- ncol(X)/M

  y.reorg <- list()
  cov.names <- list()

  for(i in 1:N){
    y.reorg[[i]] <- matrix(NA,ncol = K, nrow = M)
    cov.names[[i]] <- matrix(NA, ncol = M, nrow = M)
    for(j in 1:M){
      y.reorg[[i]][j,] <- y[[j]][i,]
      cov.names[[i]][j,] <- paste(paste(names(y)[rep(j, times = M)],i,sep = "")
                                  ,paste(names(y)[1:M],i,sep = ""), sep = ":")
    }
  }

  X.N <- list()
  Sig <- list()
  Sig.rows <- length(diag(M)[upper.tri(diag(M),diag = TRUE)])
  for(i in 1:N){
    X.N[[i]] <- X.unit[((i-1)*M+1):(i*M),]
    Sig[[i]] <- matrix(NA, nrow = Sig.rows, ncol = iters)
    Sig[[i]][,1] <- diag(M)[upper.tri(diag(M),diag = TRUE)]
  }

  ybar <- lapply(y.reorg,rowMeans)
  ybar <- do.call(rbind, ybar)
  ybar <- as.vector(t(ybar))

  Y <- do.call(rbind,y.reorg)
  Y <- as.vector(Y)

  beta <- matrix(NA, nrow = no.betas*M, ncol = iters)

  bhat <- as.matrix(solve(t(X) %*% X) %*% t(X) %*% Y)
  bhat[which(bhat < 0)] <- 1000

  beta[,1] <- bhat

  if(is.na(priors)){
  S.prior <- lapply(y.reorg,function(X){diag((apply(X,1,sd))^2)})
  nu0 <- M
  mu0 <- rep(0, times = no.betas*M)
  V0i <- diag(no.betas*M)/10000000
  }else{
    nu0 <- priors$Sig$nu0
    mu0 <- priors$beta$mu0
    V0i <- solve(priors$beta$V0)
    S.prior <- priors$Sig$S
    }

  V0imu0 <- V0i %*% mu0



  beta.temp <- rep(NA, times = no.betas*N)
  Wi.temp <- list()
  for(t in 2:iters){
    b.use <- as.vector(beta[,t-1])
    for(i in 1:N){
      xB <- as.vector(X.N[[i]] %*% b.use)
      psi <- matrix(0, ncol = M, nrow = M)
      for(k in 1:K){
        psi <- (y.reorg[[i]][,k] - xB) %*% t(y.reorg[[i]][,k] - xB) + psi
        }
      psi <- psi + S.prior[[i]]*(nu0)
      W <- rinvwishart(K + nu0, psi)
      Sig[[i]][,t] <- W[upper.tri(W, diag = TRUE)]
      Wi.temp[[i]] <- solve(W)
    }
    Oi <- as.matrix(bdiag(Wi.temp))
    XtOiX <- t(X.unit) %*% Oi %*% X.unit * K
    precis <- XtOiX + V0i
    cov.use <- solve(precis)
    bhat <- cov.use %*% (V0imu0 + t(X.unit)%*%Oi %*% ybar * K)
    beta[,t] <- as.vector(rtmvnorm(n = 1, mean = bhat[,1], sigma = cov.use, lower = rep(0, length = no.betas*M)))
    }
  Sig <- lapply(Sig,function(X, b){X[,-(1:b)]}, b = burn)
  Sig <- lapply(Sig, function(X,t){X[,seq(from = 1, to = ncol(X), by = t)]}, t = thin)

  beta <- beta[,-(1:burn)]
  beta <- beta[,seq(from = 1, to = ncol(beta), by = thin)]

  beta.return <- list()

  for(i in 1:M){
    beta.return[[i]] <- beta[((i - 1)*no.betas +1):(i*no.betas),]
  }

  names(beta.return) <- names(y)

  return(list(beta = beta.return,
              Sig = Sig,
              priors= list(beta = list(mu0 = mu0, V0 = solve(V0i)), Sig = list(nu0 = nu0, S = S.prior)), cov.structure = "location", y.cov = cov.names))
}
