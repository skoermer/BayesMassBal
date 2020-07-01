#' Two Node Process Data Simulation
#'
#' Simulates data for a stochastic two node, two component process at steady state.  Location indicies are the same as what is shown in \code{vignette("Two_Node_Process", package = "BayesMassBal")}.  Default parameters are used in an upcomming publication.
#'
#' @param k Numeric specifying the number of sample sets to be simulated.
#' @param feed List specifying qualities for the process grade.  See default for required structure.  \code{rate} is the mean feed rate.  \code{sd} is the standard deviation of the feed rate, see details.  \code{CuFeS2grade} is the mass percent CuFeS2 present in the feed.  Grade is not stochastic.
#' @param rec List specifying mean and variance in process performance.  See default for required structure. \code{rec$*component*$mean} is a vector giving mean fractional recovery of the given component for \code{c(node1,node2)}.  \code{rec$*component*$var} gives the variance in the process in a similar manner.  See details.
#' @param assayNoise List specifying standard deviations of random noise added to simulated process outputs.  See default for required structure.  The index of a vector within the list is equivalent to the index of the sampling location.  See Details section.
#' @param truncation Logical indicating if the simulation should be rerun, and previous results discgarded, until no observed values are less than 0. Default is TRUE.  See details for more information.
#'
#' @details
#'
#' Each of the \code{k} data sets collected from the \code{twonodeSim()} simulation is independent and identically distributed.
#'
#' The feed rate to the process is normally distributed with a mean of \code{feed$rate}, and a standard deviation of \code{feed$sd}.  If the feed rate is sufficiently small, and the standard deviation is sufficiently large, negative feed rates can be generated.
#'
#' Process recovery at each node is simulated from a \href{https://en.wikipedia.org/wiki/Beta_distribution}{beta distribution}, reparameterized as shown in \href{https://stats.stackexchange.com/a/12239}{this post} to make parameter specification more intiutive.  This reparameterization is only valid when \eqn{\sigma^2 \leq \mu(1-\mu)}, and the list argument \code{rec} must be specified as such.
#'
#' The simulation draws a random feed rate, draws random values for recovery, then performs process engineering calculations with these values to find mass flow rates at each sampling location for both sample components.  To these flow rates, random normally distributed noise is added, with a mean of zero and standard deviations specified in \code{assayNoise}.  If the standard deviations supplied are sufficiently large, the simulation can return negative mass flow rates.
#'
#' The argument \code{truncation = TRUE} discgards negative mass flow rates, and reruns the simulation untill all values are non-negative.  For some combinations of a large \code{k} and specifications in \code{feed} and \code{assayNoise}, this can happen frequently.  If if the simulation is run three or more times a warning will be printed that the returned expectations are unreliable.  If this is the case, expectations should be calculated using analytical or monte-carlo methods outside of the abilities of this function.  For the default parameters, can occur, but is rare.  The default parameters were chosen in a way that makes a tuncation warning highly unlikely.
#'
#' @return Returns a list of simulated data and expected values.  List members are as follows:
#' @return \item{\code{simulation}}{List of matricies giving simulated data.  \code{twonodeSim()$simulation} is structured so that it can directly be passed to the \code{\link{BMB}} function as the \code{y} argument.}
#' @return \item{\code{expectations}}{List of matricies giving expected values of the mass flow rate for each component at every location.  See the Details section for information about instances that may create reliability issues with this output.}
#'
#' @importFrom stats rnorm rbeta
#'
#' @export
#'
#' @examples
#'
#' y <- twonodeSim()$simulation
#'
#' # Then the BMB function can be run as
#' # BMB(y = y, ...)

twonodeSim <- function(k = 7, feed = list(rate = 100, sd = 6, CuFeS2grade = 1.2), rec = list(CuFeS2 = list(mean = c(98,95)/100, var = c(0.00005,0.00008)), gangue = list(mean = c(7,4)/100, var =c(0.00005,0.000025))),assayNoise = list(CuFeS2 = c(0.15,0.2,0.05,0.00005,0.005), gangue = c(5,1, 0.03, 2, 0.5)), truncation = TRUE){
  f.rate <- feed$rate
  recCu <- rec$CuFeS2$mean
  recG <- rec$gangue$mean
  g.cu <- feed$CuFeS2grade/100
  g.gangue <- 1-g.cu

  feed.mass <- rnorm(k,mean = f.rate, sd = feed$sd)
  true.cu <- true.gangue <- data.frame(matrix(NA, ncol = 5, nrow = k))
  names(true.cu) <- names(true.gangue) <- paste("B",1:5, sep = "")

  v.cu <- rec$CuFeS2$var
  v.g <- rec$CuFeS2$var

  meantemp <- c(rec$CuFeS2$mean, rec$gangue$mean)
  vartemp <- c(rec$CuFeS2$var, rec$gangue$var)

  if(!all(vartemp <= (meantemp*(1-meantemp)))){warning("Process var must be less than or equal to mean*(1-mean).")}

  rec.cu.params <- list()
  rec.cu.params[["r1"]]$alpha <- -(recCu[1]*(v.cu[1] + recCu[1]^2 - recCu[1])/v.cu[1])
  rec.cu.params[["r1"]]$beta <- (v.cu[1]^2 + recCu[1]^2 - recCu[1])*(recCu[1]-1)/v.cu[1]
  rec.cu.params[["r2"]]$alpha <- -(recCu[2]*(v.cu[2] + recCu[2]^2 - recCu[2])/v.cu[2])
  rec.cu.params[["r2"]]$beta <- (v.cu[2]^2 + recCu[2]^2 - recCu[2])*(recCu[2]-1)/v.cu[2]

  rec.g.params <- list()
  rec.g.params[["r1"]]$alpha <- -(recG[1]*(v.g[1] + recG[1]^2 - recG[1])/v.g[1])
  rec.g.params[["r1"]]$beta <- (v.g[1]^2 + recG[1]^2 - recG[1])*(recG[1]-1)/v.g[1]
  rec.g.params[["r2"]]$alpha <- -(recG[2]*(v.g[2] + recG[2]^2 - recG[2])/v.g[2])
  rec.g.params[["r2"]]$beta <- (v.g[2]^2 + recG[2]^2 - recG[2])*(recG[2]-1)/v.g[2]

  r1.cu <- rbeta(n = k, shape1 = rec.cu.params$r1$alpha, shape2 = rec.cu.params$r1$beta)
  r2.cu <- rbeta(n = k, shape1 = rec.cu.params$r2$alpha, shape2 = rec.cu.params$r2$beta)

  r1.g <- rbeta(n = k, shape1 = rec.g.params$r1$alpha, shape2 = rec.g.params$r1$beta)
  r2.g <- rbeta(n = k, shape1 = rec.g.params$r2$alpha, shape2 = rec.g.params$r2$beta)

  true.cu$B1 <- feed.mass*g.cu
  true.gangue$B1 <- feed.mass - true.cu$B1
  true.cu$B2 <- true.cu$B1 * r1.cu
  true.cu$B3 <- true.cu$B2 * r2.cu
  true.cu$B4 <- true.cu$B1 * (1-r1.cu)
  true.cu$B5 <- true.cu$B2 * (1-r2.cu)

  true.gangue$B2 <- true.gangue$B1 * r1.g
  true.gangue$B3 <- true.gangue$B2 * r2.g
  true.gangue$B4 <- true.gangue$B1 * (1-r1.g)
  true.gangue$B5 <- true.gangue$B2 * (1-r2.g)

  true.total <- true.cu + true.gangue

  s <- assayNoise$CuFeS2

  obs.cu <- t(apply(true.cu,1,noise, s = s))

  s <- assayNoise$gangue

  obs.gangue<- t(apply(true.gangue,1,noise, s = s))

  if(truncation == TRUE){
    count.while <- 0

  while(any(obs.cu < 0)  | any(obs.gangue < 0)){

    feed.mass <- rnorm(k,mean = f.rate, sd = feed$sd)

    r1.cu <- rbeta(n = k, shape1 = rec.cu.params$r1$alpha, shape2 = rec.cu.params$r1$beta)
    r2.cu <- rbeta(n = k, shape1 = rec.cu.params$r2$alpha, shape2 = rec.cu.params$r2$beta)

    r1.g <- rbeta(n = k, shape1 = rec.g.params$r1$alpha, shape2 = rec.g.params$r1$beta)
    r2.g <- rbeta(n = k, shape1 = rec.g.params$r2$alpha, shape2 = rec.g.params$r2$beta)

    true.cu$B1 <- feed.mass*g.cu
    true.gangue$B1 <- feed.mass - true.cu$B1
    true.cu$B2 <- true.cu$B1 * r1.cu
    true.cu$B3 <- true.cu$B2 * r2.cu
    true.cu$B4 <- true.cu$B1 * (1-r1.cu)
    true.cu$B5 <- true.cu$B2 * (1-r2.cu)

    true.gangue$B2 <- true.gangue$B1 * r1.g
    true.gangue$B3 <- true.gangue$B2 * r2.g
    true.gangue$B4 <- true.gangue$B1 * (1-r1.g)
    true.gangue$B5 <- true.gangue$B2 * (1-r2.g)

    true.total <- true.cu + true.gangue

    s <-  assayNoise$CuFeS2

    obs.cu <- t(apply(true.cu,1,noise, s = s))

    s <-  assayNoise$gangue

    obs.gangue<- t(apply(true.gangue,1,noise, s = s))
    count.while <- count.while + 1
  }
  if(count.while >= 2){warning("Excessive truncation has occured.  Returned expectations may be unreliable.  Reduce assayNoise values.")}
  }

  true.cu <- true.gangue <- data.frame(matrix(NA, ncol = 5, nrow = 1))
  names(true.cu) <- names(true.gangue) <- paste("B",1:5, sep = "")

  true.cu$B1 <- f.rate*g.cu
  true.gangue$B1 <- f.rate - true.cu$B1
  true.cu$B2 <- true.cu$B1 * recCu[1]
  true.cu$B3 <- true.cu$B2 * recCu[2]
  true.cu$B4 <- true.cu$B1 * (1-recCu[1])
  true.cu$B5 <- true.cu$B2 * (1-recCu[2])

  true.gangue$B2 <- true.gangue$B1 * recG[1]
  true.gangue$B3 <- true.gangue$B2 * recG[2]
  true.gangue$B4 <- true.gangue$B1 * (1-recG[1])
  true.gangue$B5 <- true.gangue$B2 * (1-recG[2])

  names(obs.cu) <-  names(obs.gangue) <- names(true.cu) <- names(true.gangue)<- paste0(rep("y", times = 5), 1:5)

  expectations <- list(CuFeS2 = t(true.cu), gangue = t(true.gangue))
  simulation <- list(CuFeS2 = t(obs.cu), gangue = t(obs.gangue))

  return(list(simulation = simulation, expectations = expectations))

}
