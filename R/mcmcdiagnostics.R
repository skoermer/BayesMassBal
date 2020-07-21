mcmcdiagnostics <- function(samps,diagnostics){

if(!is.list(diagnostics)){
  if(!is.pairlist(diagnostics)){
  diagnostics <- NULL
  }
}

SigDim <- c(length(samps$Sig),nrow(samps$Sig[[1]]))
betaDim <- c(length(samps$beta),nrow(samps$beta[[1]]))

Sigtemp <- data.frame(matrix(c(1:SigDim[2],rep(NA, times = 2*SigDim[2])), ncol = 3, nrow = SigDim[2]))
betatemp <- data.frame(matrix(c(1:betaDim[2],rep(NA, times = 2*betaDim[2])), ncol = 3, nrow = betaDim[2]))

names(Sigtemp) <- names(betatemp) <- c("index","cd", "ess")

betacd <- lapply(samps$beta, function(X, params = diagnostics){apply(X,1,function(X){do.call("geweke.diag", c(list(x = X),params))$z})})
Sigcd <- lapply(samps$Sig, function(X, params = diagnostics){apply(X,1,function(X){do.call("geweke.diag", c(list(x =X),params))$z})})

betaess <- lapply(samps$beta, function(X){apply(X,1,effectiveSize)})
Sigess <- lapply(samps$Sig, function(X){apply(X,1,effectiveSize)})

Sigreport <- list()
betareport <- list()

for(i in 1:betaDim[1]){
  name <- names(betaess)[i]
  Sigreport[[name]] <- Sigtemp
  Sigreport[[name]]$cd <- Sigcd[[i]]
  Sigreport[[name]]$ess <- Sigess[[i]]

  betareport[[name]] <- betatemp
  betareport[[name]]$cd <- betacd[[i]]
  betareport[[name]]$ess <- betaess[[i]]
}

diagnostics <- list(beta = betareport, Sig = Sigreport)
return(diagnostics)

}
