#' Point Estimate Mass Balance
#'
#' Function conducting a two component, two node, point estimate mass balance from \insertCite{willsbook}{BayesMassBal} on a two node process.  This function is provided for the purpose of comparing the performance of a Bayesian mass balance to a point estimate mass balance using any output of \code{\link{twonodeSim}}.
#' @param y A list of matrices of observed mass flow rates. Each matrix is a separate sample component.  The rows of each matrix index the sampling location, and the columns index the sample set number.  Necessary to format exactly the same as the output of \code{twonodeSim()$simulation}.  Any arguments to \code{\link{twonodeSim}} can be used.
#' @return Returns a list of vectors with mass flow rates \code{yhat} and grades \code{ahat}.  Similar format to argument \code{y}.  The index of a vector in the output is equivalent to the index of a row in \code{y}.
#'
#' @examples
#' y <- twonodeSim()$simulation
#'
#' yhat <- pointmassbal(y)$yhat
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats sd
#' @export
#'
#' @references
#' \insertRef{willsbook}{BayesMassBal}
#'
pointmassbal <- function(y){
  obs.cu <- t(y[["CuFeS2"]])
  obs.gangue <- t(y[["gangue"]])
  obs.total <- obs.cu + obs.gangue

  grade <- obs.cu/obs.total
  abar <- apply(grade,2,mean)
  var.a <- apply(grade,2,sd)^2

  Psi <- diag(var.a)

  c.init <- rep(0.5, times = 2)
  c.hat <- rep(NA, times = 2)
  V.rk1 <- var.a[1] + var.a[2]*c.init[1]^2 + (1-c.init[1])^2*var.a[4]
  V.rk2 <- var.a[2] + var.a[3]*c.init[1]^2 + (1-c.init[1])^2*var.a[5]

  c.hat[1] <- (abar[1]-abar[4])*(abar[2]-abar[4])/V.rk1
  c.hat[1] <- c.hat[1]/((abar[2]-abar[4])^2/V.rk1)

  c.hat[2] <- (abar[2]-abar[5])*(abar[3]-abar[5])/V.rk2
  c.hat[2] <- c.hat[2]/((abar[3]-abar[5])^2/V.rk2)


  tol <- 1e-6

  while(all((c.init-c.hat) > tol)){
    c.init <- c.hat

    V.rk1 <- var.a[1] + var.a[2]*c.init[1]^2 + (1-c.init[1])^2*var.a[4]
    V.rk2 <- var.a[2] + var.a[3]*c.init[1]^2 + (1-c.init[1])^2*var.a[5]

    c.hat[1] <- (abar[1]-abar[4])*(abar[2]-abar[4])/V.rk1
    c.hat[1] <- c.hat[1]/((abar[2]-abar[4])^2/V.rk1)

    c.hat[2] <- (abar[2]-abar[5])*(abar[3]-abar[5])/V.rk2
    c.hat[2] <- c.hat[2]/((abar[3]-abar[5])^2/V.rk2)

  }

  C <- matrix(c(1,-c.hat[1],0,-(1-c.hat[1]),0,0,1,-c.hat[2],0,-(1-c.hat[2])), byrow = TRUE, nrow=  2, ncol = 5)
  yield <- c.hat

  mean.feed <- mean(obs.total[,1])

  ahat <- as.vector(abar - Psi %*% t(C) %*% solve(C %*% Psi %*% t(C)) %*% C %*% abar)


  flows <- rep(mean.feed, times = 5)
  flows[2] <- flows[2]*c.hat[1]
  flows[4] <- flows[1] - flows[2]
  flows[3] <- flows[2]*c.hat[2]
  flows[5] <- flows[2]-flows[3]

  cu.flow <- ahat*flows
  gangue.flow <- (1-ahat)*flows

  ahat <- list(CuFeS2 = ahat*100, gangue = (1-ahat)*100)
  yhat <- list(CuFeS2 = cu.flow, gangue = gangue.flow)

  out <- list(yhat = yhat, ahat = ahat)
  return(out)

}
