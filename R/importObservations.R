#' Import Observed Mass Flow Rates
#'
#' Imports observed mass flow rates stored in a \code{*.csv} file and then organizes the data for use with the \code{\link{BMB}} function.
#' @param file Character string containing the name of \code{*.csv} from file which data will be read.  See Details below for valid file structures.
#' @param header Logical indicating if the first row of \code{file} file contains header information.  Current implementation of \code{importObservations} discards this information.
#' @param csv.params List of arguments to be passed to \code{\link{read.csv}}
#'
#' @details The purpose of this function is to make it easy to import and structure loosely organized data contained in a \code{*.csv} into a list for use as the \code{y} argument passed to the \code{\link{BMB}} function.
#' The entries in file must be organized as such:
#'
#' \itemize{
#' \item The first column of \code{file} must contain an integer sample location. The value of this integer must correspond to the column number used to specify linear constraints in \code{\link{constrainProcess}}.  For example, data for a given component collected at sampling location \eqn{y_2} should be indicated with a \code{2} in the first column of \code{file} used with \code{importObservations}.  In the \code{file} used with \code{\link{constrainProcess}}, the linear constraint(s) on \eqn{y_2} are indicated in the second column.
#' \item The second column of \code{file} must contain sample component names.  \strong{This field is case sensitive}.  Ensure a given sample component is named consistently, including capitalization and spacing.
#' \item Columns 3 to \eqn{K+2} of \code{file} must contain observed mass flow rates for the \eqn{K} collected sample sets.  All observations located in the same column should be collected at the same time.
#' \item Sample components of interest must be specified for each location.  If a sample component is not detected at some locations, but is detected at others, this component should be included in \code{file} with a specified mass flow rate of 0, or a very small number.
#' }
#'
#' \code{importObservations} reads the contents of \code{file}, sorts the sampling locations numerically, then creates a list of data frames.  Each data frame contains the data for a single sample component.
#'
#' @return Returns a list of data frames.  Each data frame is named according to the unique sample components specified in the second column of \code{file}.  This list object is intended to be used as the argument \code{y} for the \code{\link{BMB}} function.
#'
#' @importFrom utils read.csv write.csv
#' @importFrom stats rbeta
#'
#'
#' @export
#'
#' @examples
#'
#'  y <- importObservations(file = system.file("extdata", "twonode_example.csv",
#'                                       package = "BayesMassBal"),
#'                    header = TRUE, csv.params = list(sep = ";"))
#'
#' ## The linear constraints for this example data set are:
#' \donttest{C <- matrix(c(1,-1,0,-1,0,0,1,-1,0,-1), byrow = TRUE, ncol = 5, nrow = 2)}
#'
#' ## The X matrix for this data set can be found using:
#' \donttest{X <- constrainProcess(C = C)}

importObservations <- function(file, header = FALSE, csv.params = NULL){
  dat <- do.call(read.csv, c(list(file = file, header = header, stringsAsFactors = FALSE), csv.params))

  dat[,2] <- as.character(dat[,2])

  if(!is.character(dat[,2])){warning(paste("Second column of ", file, " must name the sample component in each row.", sep =""))}
  if(!is.integer(dat[,1])){warning(paste("First column of ", file, " must be an integer specifing the sampling location of each row.", sep = ""))}

  u.components <- unique(dat[,2])
  u.locations <- unique(dat[,1])

  y <- list()

  for(i in 1:length(u.components)){
    dat.temp <- dat[(dat[,2] == u.components[i]),]
    s <- sort(dat.temp[,1], index.return = TRUE)$ix
    dat.temp <- dat.temp[s,]

    y[[u.components[i]]] <- dat.temp[,-(1:2)]
    y[[u.components[i]]] <- as.matrix(y[[u.components[i]]])
  }

nrow.check <- lapply(y,nrow)

nrow.check <- length(unique(nrow.check))

if(nrow.check != 1){warning(paste("Number of sampling locations differ between sample components.  Check row two of ", file, " for spelling errors. See documentation for details", sep =""))}

return(y)

}
