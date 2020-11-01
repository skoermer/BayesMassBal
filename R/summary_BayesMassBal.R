#' Summary of BayesMassBal Object
#'
#' Prints a summary table containing mean values and 95% HPDI intervals for the mass flow rates, as well as the log-marginal likelihood for a \code{BayesMassBal} class object.  Options include trace plots, posterior densities, and main effects plots.
#'
#' @param object A \code{BayesMassBal} object returned from the \code{\link{BMB}} function
#' @param export  Optional character string specifying location to save a \code{*.csv} file containing summary data.  Only data related to mass flow rates is printed.
#' @param ... Additional arguments affecting the summary produced.  Not used for a \code{BayesMassBal} object.
#'
#' @details Current implementation only returns statistics for balanced mass flow rates, taken from \code{x$ybal}, and not statistics on \eqn{\beta} or variance parameters of \eqn{\sigma^2} and \eqn{\Sigma}.
#'
#' The header entry of the table \code{95\% LB} should be interpreted as the lower bound of the 95% HPDI.  Similarly, the header entry of the table \code{95\% UB} should be interpreted as the upper bound of the 95% HPDI.
#'
#' @return A summary table printed to the console, and optionally a saved \code{*.csv} file saved within the path as specified.
#'
#' @importFrom HDInterval hdi
#'
#' @export
#'

summary.BayesMassBal <- function(object, export = NA,...){

  ybal <- object$ybal
  ans <- list()

  components <- names(ybal)
  components <- c(components, "Total")
  locations <- nrow(ybal[[1]])
  ybal_total <- Reduce("+", ybal)
  ybal[["Total"]] <- ybal_total

  ans$`Mass Flow Rates` <- list()
  template_df <- data.frame(matrix(NA, ncol = 4, nrow = locations))
  names(template_df) <- c("Sampling Location", "Expected Value", "95% LB", "95% UB")
  template_df[,1] <- 1:locations

  cat("Mass Flow Rates:\n")

  for(i in 1:(length(components))){
    cat(paste("\n",components[i],":\n", sep = ""))
    temp <- template_df
    temp[,2] <- apply(ybal[[components[i]]],1,mean)
    temp[,3:4] <- unname(t(apply(ybal[[components[[i]]]],1,hdi)))
    cat(paste(c(rep("-", times = 20),"\n"), sep = "", collapse = ""))
    print(temp, row.names = FALSE)
    ans[[1]][[components[i]]] <- temp
  }

  cat("\n\nlog-marginal likelihood:\n")
  cat(paste(c(object$lml,"\n")))

  if(is.character(export)){
    csv.check <- strsplit(export, split ="[.]")[[1]]
    if(length(csv.check) == 1){
      export <- paste(export,"csv", sep = ".")
      csv.check <- strsplit(export, split ="[.]")[[1]]
    }
    if(csv.check[[2]] != "csv"){warning("Only a .csv format can be output.  Output not saved.  Check spelling of export argument", immediate. = TRUE)}
    if(length(csv.check) == 2 & csv.check[2] == ".csv"){
      export.df <- do.call("rbind",ans[[1]])
      export.df <- cbind.data.frame(`Sample Component` = rep(components, each = locations), export.df)
      write.csv(export.df, file = "test.csv", row.names = FALSE)
    }
  }else if(!is.na(export) & !is.character(export)){
    warning("\nPlease specify a character string or NA for the export argument.")
  }
}
