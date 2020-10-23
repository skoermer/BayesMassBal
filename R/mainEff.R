#' Main Effects
#'
#' Calculates the main effect of a variable, which is independent of process performance, on a function.
#'
#' @param BMBobj A \code{BayesMassBal} object originally obtained from the \code{\link{BMB}} function.  See \code{\link{BMB}}.
#' @param fn A character string naming a function with arguments of \code{BMBobj$ybal} and independent random variables \code{X}.  See Details and examples for more on function requirements.
#' @param rangex A numeric matrix.  Each column of \code{rangex} contains the minimum and maximum value of uniformly distributed random values making up vector \eqn{x}.
#' @param xj Integer indexing which element in \eqn{x} is used for conditioning for  \eqn{E_x\lbrack f(x,y)|x_j\rbrack}. If a vector is supplied the marginal main effect of each element is calculated sequentially.  The integers supplied in \code{xj} are equivalent to the indices of the columns in \code{rangex}.
#' @param N Integer specifying the length of the sequence used for \code{xj}.  Larger \code{N} trades a higher resolution of the main effect of \code{xj} for longer computation time and larger RAM requirements.
#' @param res Integer indicating the number of points to be used for each Monte-Carlo integration step.  Larger \code{res} reduces Monte-Carlo variance as the expense of computation time.
#' @param hdi.params Numeric vector of length two, used to calculate Highest Posterior Density Interval (HPDI) of the main effect \code{xj} using \code{\link[HDInterval]{hdi}}. \code{hdi.params[1] = 1} indicates \code{\link[HDInterval]{hdi}} is used, and the mean and HPDI bounds are returned instead of the every sample from the distribution of \eqn{E_x\lbrack f(x,y)|x_j\rbrack}.  The second element of \code{hdi} is passed to the \code{credMass} argument in the \code{\link[HDInterval]{hdi}} function.  The default, \code{hdi.params = c(1,0.95)}, returns 95\% HPDI bounds.
#' @param ... Extra arguments passed to the named \code{fn}
#'
#' @details
#'
#' The \code{mainEff} function returns a distribution of \eqn{E_x\lbrack f(x,y)|x_j\rbrack}, marginalized over the samples of \code{BMBobj$ybal}, giving the distribution of \eqn{E_x\lbrack f(x,y)|x_j\rbrack} which incorporates uncertainty of a chemical or particulate process.
#'
#' In the current implementation of \code{mainEff} in the \code{BayesMassBal} package, only uniformly distributed values of \eqn{x} are supported.
#'
#' The \eqn{f(x,y)} is equivalent to the supplied function named in \code{mainEff(fn)}.  For the arguments of \code{fn}, \code{ybal} is structured in a similar manner as \code{BMBobj$ybal}.  The only difference being individual columns of each matrix are used at a time, and are vectorized.  Note the way \code{ybal} is subset in the example function \code{fn_example}.  The supplied \code{X} is a matrix, with columns corresponding to each element in \eqn{x}.  The output to \code{fn} must be a vector of length \code{nrow(x)}.  The first argument of \code{fn} must be \code{X}, the second argument must be \code{BMBobj$ybal}.  Order of other arguments passed to \code{fn} through \code{...} does not matter. Look at the example closely for details!
#'
#' @return A list of \code{length(xj)} list(s).  Each list specifies output for the main effect of a \code{xj}
#' @return \item{\code{g}}{The grid used for a particular \code{xj}}
#' @return \item{\code{fn.out}}{A matrix giving results on \eqn{E_x\lbrack f(x,y)|x_j\rbrack}.  If \code{hdi.params[1] = 1}, the mean and Highest Posterior Density Interval (HPDI) bounds of \eqn{E_x\lbrack f(x,y)|x_j\rbrack} are returned.  Otherwise, samples of \eqn{E_x\lbrack f(x,y)|x_j\rbrack} are returned.  The index of each column of \code{fn.out} corresponds to the the value of \code{g} at the same index.}
#' @return \item{\code{fn}}{Character string giving the name of the function used.  Same value as argument \code{fn}.}
#' @return \item{\code{xj}}{Integer indicating the index of \eqn{x} corresponding to a grouped \code{fn.out} and \code{g}.}
#'
#' @importFrom HDInterval hdi
#'
#' @export
#'
#' @examples
#'
#' ## Importing Data, generating BMB object
#' y <- importObservations(file = system.file("extdata", "twonode_example.csv",
#'                                     package = "BayesMassBal"),
#'                    header = TRUE, csv.params = list(sep = ";"))
#'
#' C <- matrix(c(1,-1,0,-1,0,0,1,-1,0,-1), byrow = TRUE, ncol = 5, nrow = 2)
#' X <- constrainProcess(C = C)
#'
#' BMB_example <- BMB(X = X, y = y, cov.structure = "indep",
#'                    BTE = c(10,200,1), lml = FALSE, verb=0)
#'
#' fn_example <- function(X,ybal){
#'     cu.frac <- 63.546/183.5
#'     feed.mass <- ybal$CuFeS2[1] + ybal$gangue[1]
#'     ## Concentrate mass per ton feed
#'     con.mass <- (ybal$CuFeS2[3] + ybal$gangue[3])/feed.mass
#'     ## Copper mass per ton feed
#'     cu.mass <- (ybal$CuFeS2[3]*cu.frac)/feed.mass
#'     gam <- c(-1,-1/feed.mass,cu.mass,-con.mass,-cu.mass,-con.mass)
#'     f <- X %*% gam
#'     return(f)
#'     }
#'
#' rangex <- matrix(c(4.00 ,6.25,1125,1875,3880,9080,20,60,96,208,20.0,62.5),
#'                   ncol = 6, nrow = 2)
#'
#' mE_example <- mainEff(BMB_example, fn = "fn_example",rangex =  rangex,xj = 3, N = 15, res = 4)
#'
mainEff <- function(BMBobj, fn,rangex,xj,N = 50,res = 100, hdi.params = c(1,0.95),...){

  fn.name <- fn
  fn <- match.fun(fn)

  out <- list()

  if(is.null(BMBobj$ybal)){
    BMBobj$ybal <- lapply(BMBobj$beta,function(X,x){X <- x %*% X;
    row.names(X) <- paste(rep("y", times = nrow(X)), 1:nrow(X), sep ="_");
    return(X)}, x = BMBobj$X)
  }

  if (requireNamespace("tgp", quietly=TRUE)){
    LHS <- tgp::lhs
  }else{
    warning("The tgp package is required")
  }

  Ts <- ncol(BMBobj$beta[[1]])

  for(j in 1:length(xj)){

  g <- seq(from = rangex[1,xj[j]], to = rangex[2,xj[j]], length.out = N)
  exp.y <- matrix(NA, ncol = N, nrow = Ts)

  r.x <- t(rangex[,-xj[j]])

  s.v <- c((1:ncol(rangex))[-xj[j]],xj[j])
  s <- sort(s.v, index.return = TRUE)$ix

  for(i in 1:N){
  g.i <- rep(g[i], times = res)
  exp.int <- rep(NA, times = Ts)

  for(k in 1:Ts){
    U <- LHS(res,r.x)
    U <- cbind(U,g.i)[,s]
    ybal <- lapply(BMBobj$ybal,function(X,t){X[,t]},t = k)
    y <- fn(U,ybal,...)
    ey <- mean(y)
    exp.int[k] <- ey
  }
  exp.y[,i] <- exp.int
  }
  if(hdi.params[1] == 1){
   bounds <- t(apply(exp.y,2,hdi,credMass = hdi.params[2]))
   m <- apply(exp.y,2,mean)
   exp.y <- rbind(bounds[,1],m,bounds[,2])
  }

  out <- c(out,list(g = g, fn.out = exp.y,fn = fn.name, xj = xj[[j]]))

  }
  return(out)
}
