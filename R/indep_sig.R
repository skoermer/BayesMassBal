indep_sig <- function(X,y,priors,BTE = c(3000,100000,1), verb = 1, eps = sqrt(.Machine$double.eps)){

  gamdraw <- function(x, shape = alpha){
    #rgamma(1,shape = shape, rate = x)
    rinvgamma(1,shape = shape, scale = x)
  }
  nDigits <- function(x){
    truncX <- floor(abs(x))
    if(truncX != 0){
      floor(log10(truncX)) + 1
    } else {
      1
    }
  }

  # Number of tests
  K <- ncol(y[[1]])
  ## Sample Locations
  N <- nrow(y[[1]])
  ## Components
  M <- length(y)

  burn <- BTE[1]
  iters <- BTE[2]
  thin <- BTE[3]

  no.betas <- ncol(X)

  x.unit <- X
  X.unit <- as.matrix(bdiag(replicate(M,x.unit,simplify = FALSE)))
  X <- do.call("rbind", replicate(K, X.unit, simplify=FALSE))

  X.M <- list()
  Sig <- list()

  for(i in 1:M){
    X.M[[i]] <- X.unit[((i-1)*N+1):(i*N),]
    Sig[[i]] <- matrix(NA, nrow = N, ncol = iters)
  }

  ybar <- lapply(y,rowMeans)

  Y <- do.call(cbind,y)
  Y <- as.vector(Y)

  beta <- replicate(M, matrix(NA, nrow = no.betas, ncol = iters), simplify = FALSE)

  names(beta) <- names(Sig) <- names(y)

  bhat <- as.vector(solve(t(X) %*% X) %*% t(X) %*% Y)

  phi.priors <- list()
  bhat[which(bhat < 0)] <- 0

  if(is.na(priors)){
    b0 <- 0.000001
    a0 <- 0.000001
    mu0 <- bhat
    dgts <- sapply(bhat,nDigits)
    V0i <- diag(1/(10^(dgts+6)))
    for(i in 1:M){
      phi.priors[[i]] <- cbind.data.frame(a = rep(a0, times = N), b = rep(b0, times = N))
    }
  }else{
    a0 <- priors$phi$a
    b0 <- priors$phi$b
    mu0 <- priors$beta$mu0
    V0i <- solve(priors$beta$V0)
  }

  bhat[which(bhat == 0)] <- eps
  b.use <- bhat

  V0imu0 <- V0i %*% mu0


  alpha <- (K)/2 + a0
  rate <- matrix(NA, ncol= M, nrow =N)

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
      for(j in 1:N){
        rate[j,i] <- 0.5*(t(xB[j] - y[[i]][j,]) %*% (xB[j] - y[[i]][j,]) )+ b0
        }
      Sig[[i]][,t] <- s.temp <- sapply(rate[,i], gamdraw, shape = alpha)
      Wi <- diag(1/s.temp)
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

  samps <- list(beta = beta,
                Sig = Sig,
                priors= list(beta = list(mu0 = mu0, V0 = solve(V0i)), phi = phi.priors), cov.structure = "indep")

  return(samps)
}
