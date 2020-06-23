#' Mass Balance Data
#'
#' Imports observed mass flow rates stored in a \code{*.csv} file and then organizes the data for use with the \code{\link{bayes.massbal}} function.
#' @param file Character string containing the name of \code{*.csv} from file which data will be read.  See Details below for valid file structures.
#' @param header Logical indicating if the first row of \code{file} file contains header information.  Current implementation of \code{massbal.data} discards this information.
#'
#' @details The purpose of this function is to make it easy to import and structure loosly organized data contained in a \code{*.csv} into a list for use with \code{\link{bayes.massbal}}.
#' The entries in file must be organized as such:
#'
#' \itemize{
#' \item The first column of \code{file} must contain an integer sample location which corresponds to the column number used for linear constraints.  For example, data for a given component collected at sampling location \eqn{y_2} should be indicated with a \code{2} in the first column of \code{file} used with \code{massbal.data}.  In the \code{file} used with \code{\link{constrain.process}}, the linear constraint on \eqn{y_2} is indicated in the second column.
#' \item The second column of \code{file} must contain sample component names.  \strong{This field is case sensitive}.  Ensure a given sample component is named consistently, including capitilazation and spacing.
#' \item Columns 3 to \eqn{K+2} of \code{file} must contain observed mass flow rates for the \eqn{K} mass flow rate observations.  All observations located in the same column should be collected at the same time.
#' \item Sample components of interest must be specified for each location.  If a sample component is not detected at some locations, but is detected at others, this component should be included in \code{file} with a specified mass flow rate of 0, or a very small number.
#' }
#'
#' \code{massbal.data} reads the contents of \code{file}, sorts the sampling locations numerically, then creates a list of data frames.  Each data frame contains the data for a single sample component.
#'
#' @return Returns a list of data frames.  Each data frame is named according to the unique sample components specified in the second column of \code{file}.  This list object is intended to be used with the \code{\link{bayes.massbal}} function.
#'
#' @importFrom utils read.csv write.csv
#' @importFrom stats rbeta
#'
#'
#' @export
#'
#' @examples
#'
#' # Data is simulated for a single node,
#' # one input, two output, two component process,
#' # with three sets of observations.
#'
#' y.save <- data.frame(matrix(NA, ncol = 5, nrow=  6))
#'
#' # Specify the sampling location
#' y.save[,1] <- rep(1:3, times = 2)
#'
#' # Specify the sample component
#'
#' y.save[,2] <- rep(c("CuFeS2","gangue"), each = 3)
#'
#' # Then the data is simulated, shuffled, and saved to be used as if it is real user collected data
#'
#' f.rate <- 100
#' f.rate.cu <- 1.2/100 * f.rate
#' f.rate.gangue <- (100-1.2)/100 *f.rate
#' cu.rec <- rbeta(3,shape1 = 5, shape2 = 2)
#' gangue.rec <- rbeta(3, shape1 = 2, shape2 =5)
#'
#' y.save[1,3:5] <- f.rate.cu
#' y.save[2,3:5] <- cu.rec * f.rate.cu
#' y.save[3,3:5] <- (1-cu.rec) * f.rate.cu
#'
#' y.save[4, 3:5] <- f.rate.gangue
#' y.save[5,3:5] <- f.rate.gangue * gangue.rec
#' y.save[6, 3:5] <- f.rate.gangue * (1-gangue.rec)
#'
#' # shuffeling the generated data
#'
#' y.save <- y.save[sample(1:nrow(y.save), size = nrow(y.save)),]
#'
#' # writing the .csv file
#'
#' write.csv(y.save, file = "example_file.csv", col.names = FALSE, row.names = FALSE)
#'
#' # Using the massbal.data function to read the data and
#' # organize it into the list required for bayes.massbal
#'
#' y <- massbal.data(file = "example_file.csv")

massbal.data <- function(file, header = FALSE){
  dat <- read.csv(file, header =header, stringsAsFactors = FALSE)

  if(!is.character(dat[,2])){warning(paste("Second column of ", file, " must name the sample component in each row.", sep =""))}
  if(!is.integer(dat[1,])){warning(paste("First column of ", file, " must be an integer specifing the sampling location of each row.", sep = ""))}

  u.components <- unique(dat[,2])
  u.locations <- unique(dat[,1])

  y <- list()

  for(i in 1:length(u.components)){
    dat.temp <- dat[(dat[,2] == u.components[i]),]
    s <- sort(dat.temp[,1], index.return = TRUE)$ix
    dat.temp <- dat.temp[s,]

    y[[u.components[i]]] <- dat.temp
  }

nrow.check <- lapply(y,nrow)

nrow.check <- length(unique(nrow.check))

if(nrow.check != 1){warning(paste("Number of sampling locations differ between sample components.  Check row two of ", file, " for spelling errors. See documentation for details", sep =""))}

return(y)

}
