scale <- function(y,X, cov.structure){

  Xold <- X
  K <- ncol(y[[1]])
  N <- nrow(y[[1]])
  M <- length(y)

  if(cov.structure == "location"){

    X.unit <- bdiag(replicate(M,X,simplify = FALSE))
    X.unit <- as.matrix(X.unit)
    X.shuffle <- matrix(1:(N*M), ncol = N, nrow = M, byrow = TRUE)
    X.shuffle <- as.vector(X.shuffle)
    X.unit <- X.unit[X.shuffle,]
    y.reorg <- list()
    X.N <- list()
    for(i in 1:N){
      X.N[[i]] <- X.unit[((i-1)*M+1):(i*M),]
      y.reorg[[i]] <- matrix(NA,ncol = K, nrow = M)
      for(j in 1:M){
        y.reorg[[i]][j,] <- y[[j]][i,]
      }
    }
    y.use <- y.reorg
  }else{y.use <- y
  X.N <- replicate(M,Xold, simplify = FALSE)
  }

  ymeans <- lapply(y.use, function(X){rowMeans(X)})
  redfctr <- lapply(ymeans, function(X){ftemp <- signif(X/1000, digits = 1);
  ftemp[ftemp <= 1] <- 1;
  return(ftemp)})
  incfctr <- lapply(ymeans, function(X){ftemp <- signif(X/0.001, digits = 1);
  ftemp[ftemp >= 1] <- 1
  ftemp[ftemp == 0] <- 1;
  return(ftemp)})

  A <- list()
  for(i in length(y.use)){
    y.use[[i]] <- y.use[[i]]/redfctr[[i]]
    y.use[[i]] <- y.use[[i]]/incfctr[[i]]

    X.N[[i]] <- t(t(X.N[[i]])/redfctr[[i]])
    X.N[[i]] <- t(t(X.N[[i]])/incfctr[[i]])

    Atemp <- diag(nrow(X.N[[i]]))
    Atemp <- Atemp*redfctr[[i]]
    Atemp <- Atemp*incfctr[[i]]
    A[[i]] <- Atemp
  }


  return(list(yscale = y.use, Xscale = X.N, A = A))

}
