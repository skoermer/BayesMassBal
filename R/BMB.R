#' Bayesian Mass Balance
#'
#' Allows the user to specify the covariance structure for a Bayesian mass balance, simulates draws from reconciled masses and relevant covariance matrix, and approximates the log-marginal likelihood.
#' @param X A matrix that maps constrained masses to observed masses.  Can be built from the function \code{\link{constrainProcess}}, see documentation for details.
#' @param y A list of matrices of observed mass flow rates. Each matrix is a separate sample component.  The rows of each matrix index the sampling location, and the columns index the sample set number.  Can be specified using the \code{\link{importObservations}} function.
#' @param cov.structure Character string. \code{"indep"} allows for no correlation.  \code{"component"} indicates correlation within an individual sample component.  \code{"location"} indicates correlation within an individual sampling location.  Not specifying \code{cov.structure} defaults to the \code{"indep"} structure.
#' @param priors List or character string. When the default value \code{priors = "default"} is used, the \code{BMB} uses a set of default conjugate priors.  When passing a list to the argument, the list must contain user specified hyperparameter values for each conjugate prior.  To see the required list structure run \code{BMB} with \code{BTE = c(1,2,1)} and inspect the output.  When \code{priors = "Jeffreys"} the Jeffreys priors for \eqn{\Sigma} and \eqn{\sigma^2} with known mean given in \insertCite{priorlist}{BayesMassBal}.  When \code{priors = "Jeffreys"}, the prior used for \eqn{\beta} is proportional to the indicator function \eqn{I\lbrack \beta > 0 \rbrack}.  See Details for more information.
#' @param BTE Numeric vector giving \code{c(Burn-in, Total-iterations, and Every)} for MCMC approximation of target distributions. The function \code{BMB} produces a total number of samples of \eqn{(T - B)/E}.  \eqn{E} specifies that only one of every \eqn{E} draws are saved. \eqn{E > 1} reduces autocorrelation between obtained samples at the expense of computation time.
#' @param lml Logical indicating if the log-marginal likelihood should be approximated.  Default is \code{FALSE}, which reduces computation time.  Log-marginal likelihood is approximated using methods in \insertCite{chib}{BayesMassBal}.
#' @param ybal Logical indicating if the mass balanced samples for each \eqn{y} should be returned.  Default is \code{TRUE}. Setting \code{ybal=FALSE} results in a savings in RAM and computation time.
#' @param diagnostics Logical or list indicating if diagnostic functions \code{\link[coda]{geweke.diag}} and \code{\link[coda]{effectiveSize}} \insertCite{coda}{BayesMassBal} should computed for the obtained samples.  The default of \code{TRUE} indicates diagnostics should be run with their default parameters.  Alternatively, passing a list of the structure \code{list(frac1 = 0.1, frac2 = 0.5)} will run both diagnostics and allow \code{\link[coda]{geweke.diag}} to be run with parameters other than the default.
#' @param verb Numeric indicating verbosity of progress printed to R-console.  The default of 1 prints messages and a progress bar to the console during all iterative methods. \code{verb = 0} indicates no messages are printed.
#' @return Returns a list of outputs
#' @return \item{\code{beta}}{List of matrices of samples from the distribution of reconciled data.  Each matrix in the list is a separate sample component. Each column of a matrix in  \code{beta} is a draw from the target distribution.}
#' @return \item{\code{Sig}}{List of matrices containing draws from the distribution of each covariance matrix.  If \code{S.t} is the \eqn{t^{th}} draw from the distribution of covariance matrix \code{S} and:
#' \itemize{\item \code{cov.structure = "indep"}, the \eqn{t^{th}} column of a matrix in \code{Sig} is \code{diag(S.t)}.
#' \item\code{cov.structure = "component"} or \code{"location"}, the \eqn{t^{th}} column of a matrix in \code{Sig} is equal to \code{S.t[upper.tri(S.t, diag = TRUE)]}.}}
#' @return \item{\code{priors}}{List of prior hyperparameters used in generating conditional posterior distributions and approximating log-marginal likelihood.  The structure of the input argument \code{priors} is required to be the same as the structure of this returned list slice.  See Details.}
#' @return \item{\code{cov.structure}}{Character string containing the covariance structure used.}
#' @return \item{\code{y.cov}}{List of character matrices indicating details for the structure of each covariance matrix.  Only returned when \code{cov.structure = "location"}}
#' @return \item{\code{lml}}{Numeric of the log-marginal likelihood approximation. Returns \code{NA} when \code{lml = FALSE}}
#' @return \item{\code{diagnostics}}{List containing results from diagnostic functions \code{\link[coda]{geweke.diag}} and \code{\link[coda]{effectiveSize}}}
#' @return \item{\code{ybal}}{List of samples from the distribution of reconciled mass flow rates, in the same format as the function argument \code{y}.  Produced with argument \code{ybal = TRUE}.  Equivalent to \code{lapply(BMB(...)$beta,function(X,x){x \%*\% X} , x = X)} .  Viewing this output is more intuitive than viewing samples of \code{beta}, at the expense of RAM and some computation time.}
#' @return \item{\code{X}}{The function argument \code{X} is passed to the output so that it can be used with other \code{BayesMassBal} functions.}
#' @return \item{\code{type}}{Character string used by \code{\link{plot.BayesMassBal}}.  \code{type = "BMB"} for an object returned from the \code{BMB} function.}
#'
#' @details
#'
#' See \code{vignette("Two_Node_Process", package = "BayesMassBal")} for further details on how to use function outputs.
#'
#' When the \code{priors} argument is left unspecified, a set of default conjugate priors are used, which are chosen to allow \code{BMB()} to work well in a general setting.  In the current version of the \code{BayesMassBal} package, only the conjugate priors stated below can be used, but hyperparameter values can be specified by the user using the \code{priors} argument.
#'
#' The prior distribution on \code{beta} is a normal distribution truncated at 0.  The mean of this distribution before truncation is the \href{https://en.wikipedia.org/wiki/Ordinary_least_squares}{ordinary least squares} (OLS) estimate of \eqn{\beta}. OLS estimates less than 0, are changed to 0.  The prior variance, before truncation, of each element of \eqn{\beta} is set to:
#'
#' \deqn{10^{\mathrm{number of integer digits of an element of } \beta + 6}}
#'
#' Currently, there is only support for diagonal prior covariance matrices for \eqn{\beta}
#'
#' When \code{cov.structure = "indep"} the error of all observations in a sample set are independent. An \href{https://en.wikipedia.org/wiki/Inverse-gamma_distribution}{inverse gamma} prior distribution, with \eqn{\alpha_0 = 0.000001} and \eqn{\beta_0 = 0.000001}, is placed on the variance of the mass flow rate for each sample component at each sample location.
#'
#' When \code{cov.structure = "component"} or \code{"location"}, the prior distribution on \eqn{\Sigma_i} is \href{https://en.wikipedia.org/wiki/Inverse-Wishart_distribution}{inverse Wishart} \eqn{(\nu_0, \nu_0 \times S_0)}. The degrees of freedom parameter, \eqn{\nu_0}, is equal to the dimension of \eqn{\Sigma_i}.  The scale matrix parameter is equal to a matrix, \eqn{S_0}, with the sample variance of the relevant observation on the diagonal, multiplied by \eqn{\nu_0}.
#'
#' The user is able to specify the prior hyperparameters of the mean and variance of \code{beta}, \eqn{\alpha_0} and \eqn{\beta_0} for each \eqn{\sigma^2}, and the degrees of freedom and scale matrix for each \eqn{\Sigma_i} using the \code{priors} argument.  It is advisable for the user to specify their own prior hyperparameters for \eqn{p(\sigma^2)} if the variance of any element is well under 1, or \eqn{p(\beta)} if the there is a wide range in the magnitude of observations.
#'
#' When \code{priors = "Jeffreys"} \href{https://en.wikipedia.org/wiki/Jeffreys_prior}{Jeffreys} priors are used for the prior distribution of the variance and covariance parameters.  Priors used are \eqn{p(\sigma^2) \propto \frac{1}{\sigma^2}} and \eqn{p(\Sigma) \propto |\Sigma|^{-(p+1)/2}}, as listed in \insertCite{priorlist}{BayesMassBal}.  The Jeffreys prior for a \eqn{\beta} with infinite support is \eqn{p(\beta) \propto 1}.  To preserve the prior information that \eqn{\beta > 0}, \eqn{p(\beta)\propto I\lbrack \beta > 0 \rbrack} is chosen.  It is not possible to calculate log-marginal likelihood using the methods in \insertCite{chib}{BayesMassBal} with Jeffreys priors.  Therefore, if \code{priors = "Jeffreys"} and \code{lml = TRUE}, the \code{lml} argument will be ignored and a warning will be printed.
#'
#' \code{lml} is reported in base \eqn{e}.  See \href{https://en.wikipedia.org/wiki/Bayes_factor#Interpretation}{here} for some guidance on how to interpret Bayes Factors, but note log base 10 is used on Wikipedia.
#'
#' @examples
#' y <- importObservations(file = system.file("extdata", "twonode_example.csv",
#'                                     package = "BayesMassBal"),
#'                    header = TRUE, csv.params = list(sep = ";"))
#'
#' C <- matrix(c(1,-1,0,-1,0,0,1,-1,0,-1), byrow = TRUE, ncol = 5, nrow = 2)
#' X <- constrainProcess(C = C)
#'
#' BMB_example <- BMB(X = X, y = y, cov.structure = "indep",
#'                    BTE = c(10,300,1), lml = FALSE, verb=0)
#'
#' summary(BMB_example)
#'
#' @importFrom Rdpack reprompt
#' @importFrom Matrix bdiag
#' @importFrom tmvtnorm rtmvnorm dtmvnorm
#' @importFrom LaplacesDemon rinvwishart dinvwishart rinvgamma dinvgamma
#' @importFrom stats sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom coda geweke.diag effectiveSize
#' @export
#'
#' @references
#' \insertRef{chib}{BayesMassBal}
#' \insertRef{gibbsexpl}{BayesMassBal}
#' \insertRef{coda}{BayesMassBal}
#' \insertRef{priorlist}{BayesMassBal}

BMB <- function(X, y, cov.structure = c("indep","component","location"), priors = "default", BTE = c(500,20000,1), lml = FALSE, ybal = TRUE, diagnostics = TRUE, verb = 1){

  if(all(cov.structure == c("indep","component","location"))){cov.structure <- "indep"}
  if(is.null(names(y))){names(y) <- paste0("component",1:length(y))}
  if(!is.matrix(y[[1]])){y <- lapply(y,as.matrix)}
  if(priors == "Jeffreys" & lml == TRUE){
    lml <- FALSE
    warning("Methods used do not allow lml approximation when argument: priors = \"Jeffreys\".  lml has been set to FALSE\n", immediate. = TRUE)
  }

  if(cov.structure == "indep"){
    samps <- indep_sig(X = X,y = y,priors = priors,BTE = BTE, verb = verb)
    chib.out <- NA
    if(lml == TRUE){chib.out <- chib_indep(s.l = samps, X = X, y = y, verb = verb)$lml}
  }else if(cov.structure == "component"){
    samps <- component_sig(X =X,y = y,priors = priors,BTE = BTE, verb = verb)
    chib.out <- NA
    if(lml == TRUE){chib.out <- chib_component(s.l = samps,X = X, y = y, verb = verb)$lml}
  }else if(cov.structure == "location"){
    samps <- location_sig(X = X, y = y, priors = priors, BTE = BTE, verb = verb)
    chib.out <- NA
    if(lml == TRUE){chib.out <- chib_location(s.l = samps, X = X, y = y,verb = verb)$lml}
  }else{warning("Please select a valid covariance structure.  See Argument cov.structure in documentation.", immediate. = TRUE)}

  samps$lml <- chib.out

  if(any(diagnostics != FALSE)){
    diagnostics <- mcmcdiagnostics(samps,diagnostics)
    samps$diagnostics <- diagnostics
  }

  if(ybal == TRUE){
    samps$ybal <- lapply(samps$beta,function(X,x){X <- x %*% X;
    row.names(X) <- paste(rep("y", times = nrow(X)), 1:nrow(X), sep ="_");
    return(X)}, x = X)
  }
  samps$X <- X
  samps$type <- "BMB"
  class(samps) <- "BayesMassBal"
  if(verb != 0){message("Done!")}
  return(samps)
}

