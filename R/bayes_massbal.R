#' Bayesian Mass Balance
#'
#' One function that allows the user to specify covariance structure for a Bayesian mass balance, simulates draws from reconsciled masses and relevant covariance matrix, and approximates the log(marginal likelihood).
#' @param X A matrix that maps constrained masses to observed masses.  Can be built from the function \code{\link{constrain.process}}, see documentation for details.
#' @param y A list of matricies of observed mass flow rates. Each matrix is a seperate sample component.  The rows of each matrix index sampling location, and the columns index sample set number.
#' @param cov.structure Character string. \code{"indep"} allows for independent error.  \code{"component"} indicates error for a sample component is correlated, and error between components is independent.  \code{"location"} indicates correlated error between all sample components at a given location, with independent error between location.  Not specifying \code{cov.structure} defaults to the \code{"indep"} structure.
#' @param priors Optional list of user specified hyperparameters for conjugate priors.  To see the required list structure run \code{bayes.massbal} with \code{BIT = c(1,2,1)} and inspect the output.  When not specified, the function defaults to the priors: \eqn{p(\beta) = } Truncated Normal \eqn{(0,10000000)}, for indpendent error \eqn{p(\phi_{i,j}) = } Gamma \eqn{(1,50000000)}, for location and component error \eqn{p(\Sigma_i)} Inverse Wishart\eqn{(\nu_0,\nu_0 \times S_0)} where \eqn{S_0} is a diagonal matrix containing the emperical variance of the observations for each location and component, and \eqn{\nu_0} is equal to the dimension of \eqn{\Sigma_i}.
#' @param BIT Numeric vector giving \code{c(Burn-in, total-Iterations, and Thinning)} for MCMC approximation of target distributions. The function \code{bayes.massbal} produces a total number of samples of \eqn{(I - B)/T}.  Thinning reduces autocorrelation between consecutive samples at the expense of computation time.
#' @param lml Logical indicating if the log marginal likelihood should be approximated.  Default is \code{FALSE}, which reduces computation time.  Log-marginal likekihood is approximated using methods in \insertCite{chib}{BayesMassBal}.
#' @return Returns a list of outputs
#' @return \item{\code{beta}}{List of matricies of samples from the distribution of reconsciled data.  Each matrix in the list is a seperate sample component.}
#' @return \item{\code{Sig}}{List of matricies containing draws from the distribution of each covariance matrix.  If \code{cov.structure = "indep"} the elements of each matrix are the diagonal elements of the covariance between sample locations for a given component.  For \code{cov.structure = "component"} or \code{"location"}, if \code{S.t} is a draw from the distribution of covariance matrix \code{S}, the \eqn{t^th} column in a listed matrix is equal to \code{S.t[upper.tri(W, diag = TRUE)]}}
#' @return \item{\code{priors}}{List of prior hyperparameters used in generating conditional posterior distributions and approximating log marginal likelihood.  The structure of the input variable \code{priors} is required to be the same as the structure of the returned prior hyperparameters.}
#' @return \item{\code{cov.structure}}{Character string containing the covariance structure used.}
#' @return \item{\code{y.cov}}{List of character matricies indicating details for the structure of each covariance matrix.  Only returned when \code{cov.structure = "location"}}
#' @return \item{\code{lml}}{Numeric of the log marginal likelihood approximation. Returns \code{NA} when \code{lml = FALSE}}
#'
#' @importFrom Rdpack reprompt
#' @importFrom Matrix bdiag
#' @importFrom tmvtnorm rtmvnorm dtmvnorm
#' @importFrom LaplacesDemon rinvwishart dinvwishart
#' @importFrom stats dgamma rgamma sd
#' @export
#'
#' @references
#' \insertRef{chib}{BayesMassBal}
#' \insertRef{gibbsexpl}{BayesMassBal}

bayes.massbal <- function(X,y,cov.structure = c("indep","component","location"),
                          priors = NA,BIT = c(500,20000,1), lml = FALSE){

  if(cov.structure == "indep" | all(cov.structure == c("indep","component","location"))){
    samps <- indep.sig(X = X,y = y,priors = priors,BIT = BIT)
    chib.out <- NA
    if(lml == TRUE){chib.out <- chib.indep(s.l = samps, X = X, y = y)$lml}
  }else if(cov.structure == "component"){
    samps <- component.sig(X =X,y = y,priors = priors,BIT = BIT)
    chib.out <- NA
    if(lml == TRUE){chib.out <- chib.component(s.l = samps,X = X, y = y)$lml}
  }else if(cov.structure == "location"){
    samps <- location.sig(X = X, y = y, priors = priors, BIT = BIT)
    chib.out <- NA
    if(lml == TRUE){chib.out <- chib.location(s.l = samps, X = X, y = y)$lml}
  }else{warning("Please select a valid covariance structure.  See variable cov.structure in documentation.", immediate. = TRUE)}
  samps$lml <- chib.out
  return(samps)
}

