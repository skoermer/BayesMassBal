#' Matrix Constraining Process
#'
#' Generates matrix \eqn{X} which maps constrained masses \eqn{\beta} to observed masses \eqn{y} for an individual sample component, when given a linear system of constraining equations.
#' @param C Matrix of constraints for a process at steady state.    See Details below.
#' @param file Character string indicating file path for a \code{*.csv} file containing linear constraints.  Only values of -1, 0, and 1 are valid.  The first row in the \code{file} is required to be a header naming the sampling locations.
#' @details The output of this function is meant to be used as the input parameter \code{X} in \code{\link{BMB}}.  The matrix \code{C}, or imported matrix from \code{file}, indexes sampling locations via columns, and number of constraints via rows. Only values of -1, 0, and 1 are valid, and indicate mass leaving a node, a location that is not relevant to a node, and mass entering a node respectively.  Constraints only need to be indicated around each node.  Additional constraints are redundant.  Each sample component is subject to the same constraints, and therefore the constraints given to \code{constrain.process} do not need to be repeated for each component.  A quick example, \code{C} for a single process, where the constraint is \eqn{y_1 = y_2 + y_3}, and there are two sample components, would be \code{matrix(c(1,-1,-1), nrow=  1,ncol = 3}).
#' @return Returns the matrix \eqn{X} which maps \eqn{\beta} to observed masses \eqn{y}.  No changes need to be made to \code{X} when using with \code{\link{BMB}}.
#' @importFrom pracma rref
#' @importFrom utils read.csv
#' @export
#' @examples
#' # For a 2 node process with 1 input, 3 outputs
#'
#' C <- matrix(c(1,-1,0,-1,0,0,1,-1,0,-1), byrow = TRUE, ncol = 5, nrow = 2)
#' constrain.process(C = C)
constrain.process <- function(C = NULL, file = FALSE){
  if(is.character(file)){
    C <- as.matrix(read.csv(file))
  }
  if(!all(C %in% c(-1,0,1))){
    warning("Constraint matrix must contain only the values -1, 0, or 1, indicating negative, null, or positive mass flows.")
  }


  Cc <- rref(C)
  dim.beta <- sum(colSums(Cc) < 0)
  indep.loc <- which(colSums(Cc) < 0)
  dep.loc <- which(colSums(Cc) > 0)

  if((length(dep.loc)+length(indep.loc)) != ncol(C)){
    warning("Constraint matrix is not valid.  See package documentation.", immediate. = TRUE)
  }

  X <- matrix(NA, ncol = length(indep.loc), nrow = ncol(C))
  indep.mat <- diag(length(indep.loc))
  X[indep.loc,] <- indep.mat
  X[dep.loc,] <- -Cc[,indep.loc]

  return(X)
}
