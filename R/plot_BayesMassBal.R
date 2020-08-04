#' Plots BayesMassBal Object
#'
#' Visualizes data from a \code{BayesMassBal} class object in a user specified way.  Options include trace plots, posterior densities, and main effects plots.  Meant to be a quick diagnostic tool, and not to produce publication quality plots.
#'
#' @param x A \code{BayesMassBal} object returned from the \code{\link{BMB}} function
#' @param sample.params List to be used for indicating model parameter samples used for creation of plot(s).  See details for required structure.
#' @param layout Character string indicating the desired data to be plotted.  \code{"trace"} produces trace plots of sequential parameter draws. \code{"dens"} produces densities of posterior draws.
#' @param hdi.params Numeric vector of length two, used to draw Highest Posterior Density Intervals (HPDI) using \code{\link[HDInterval]{hdi}}, and otherwise ignored. \code{hdi.params[1] = 1} indicates \code{\link[HDInterval]{hdi}} bounds should be drawn.  The second element of \code{hdi} is passed to \code{credMass} in the \code{\link[HDInterval]{hdi}} function.  The default, \code{hdi.params = c(1,0.95)}, plots the 95\% HPDI bounds.
#' @param ... Passes extra arguments to \code{plot()}
#'
#' @details
#'
#' The list of \code{sample.params} requires a specific structure dependent on the choice of \code{layout} and the desired plots.
#'
#' If \code{layout = "trace"} or \code{layout = "dens"}, \code{names(list)} must contain each model parameter desired for plotting.  The structure under the model parameter names must be the same as to the structure of the relevant subset of the \code{BayesMassBal} object to be used.  For example, if a \code{BayesMassBal} object is created using a process with sample components \code{c("CuFeS2","gangue")} and the users wants plots of reconciled masses \eqn{y_1} and \eqn{y_2} for both components to be created, \code{params = list(y.bal = list(CuFeS2 = c(1,2), gangue = c(1,2))} should be used.  Note, \code{str(params)} mimics \code{str(x)}, while the vectors listed simply index the desired model parameters to be plotted.
#'
#' See \code{vignette("Two_Node_Process", package = "BayesMassBal")} for an example of the required structure.
#'
#' @return Plots \code{BayesMassBal} object based on arguments passed to \code{plot}.
#'
#' @importFrom HDInterval hdi
#' @importFrom graphics plot par abline plot.new
#' @importFrom stats density
#'
#' @export


plot.BayesMassBal <- function(x,sample.params = NA,layout = c("trace","dens"),hdi.params = c(1,0.95),...){

  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))

  sample.names <- names(sample.params)
  samples <- list()
  plot.names <- character()

  for(i in 1:length(sample.names)){
    sample.subset.names <- names(sample.params[[sample.names[i]]])
    sample.subset <- list()

    for(j in 1:length(sample.subset.names)){
      object.use <- sample.params[[sample.names[i]]][[sample.subset.names[j]]]
      plot.names <- c(plot.names,paste(sample.names[i], object.use, sample.subset.names[j], sep = "_"))
      sample.subset[[j]] <- x[[sample.names[i]]][[sample.subset.names[j]]][object.use ,]
    }
    ## make sure to keep names
    samples[[sample.names[i]]] <- do.call(rbind, sample.subset)
  }

  samples <- do.call(rbind,samples)

  nplots <- nrow(samples)

  nrow.layout <- floor(sqrt(nplots))
  ncol.layout <- ceiling(nplots/nrow.layout)
  plot.spaces <- nrow.layout * ncol.layout

  if(hdi.params[1] == 1){
  hdpi <- apply(samples,1,function(X,pct = hdi.params[2]){hdi(X,pct)})
  }else{hdpi <- matrix(NA, ncol = 2, nrow = nplots)}


  if(layout == "trace"){
    layout(mat = matrix(1:plot.spaces,nrow = nrow.layout, ncol = ncol.layout, byrow = TRUE))
    par(mar = c(4,4,2,1))
    for(i in 1:nplots){
      plot(samples[i,], type = "l", ylab = plot.names[i], ...)
      abline(h = hdpi[,i], col = "red")
      abline(h = mean(samples[i,]), col = "darkgreen")
    }
    if((plot.spaces - nplots) > 0){
    for(i in 1:(plot.spaces-nplots)){
      plot.new()
    }
   }
  }else if(layout == "dens"){
    d.use <- apply(samples,1,function(X){density(X, from = mean(X)-3.5*sd(X), to = mean(X) + 3.5*sd(X))})
    par(mar = c(4,4,2,1))
    for(i in 1:nplots){
      plot(d.use[[i]],xlab = plot.names[i], ylab = "", main = "", ...)
      abline(v = hdpi[,i], col = "red")
      abline(v = mean(samples[i,], col = "darkgreen"))
    }
    if((plot.spaces - nplots) > 0){
      for(i in 1:(plot.spaces-nplots)){
        plot.new()
      }
    }

  }



}
