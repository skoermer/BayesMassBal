#' Plots BayesMassBal Object
#'
#' Visualizes data from a \code{BayesMassBal} class object in a user specified way.  Options include trace plots, posterior densities, and main effects plots.  Meant to be a quick diagnostic tool, and not to produce publication quality plots.
#'
#' @param x A \code{BayesMassBal} object returned from the \code{\link{BMB}} function
#' @param sample.params List to be used for indicating model parameter samples used for creation of plot(s).  See details for required structure.
#' @param layout Character string indicating the desired data to be plotted.  \code{"trace"} produces trace plots of sequential parameter draws. \code{"dens"} produces densities of posterior draws.  Argument ignored when \code{x$type = "time-series"}.
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
#' @importFrom graphics plot par abline plot.new legend text
#' @importFrom stats density
#'
#' @export


plot.BayesMassBal <- function(x,sample.params = NA,layout = c("trace","dens"),hdi.params = c(1,0.95),...){

  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))

  if(x$type == "time-series"){
    c <- c("#E87722","#75787B")
    l.wid <- 1
    y <- x$y
   if(!is.null(x$samples$expectation)){

     mean.exp <- mean(x$samples$expectation)
     mean.alpha <- mean(x$samples$alpha)

     r.alpha<- range(x$samples$alpha)
     d.alpha <- density(x$samples$alpha, from = r.alpha[1], to = r.alpha[2])



     leg.lab <- c("Data","Expected Steady State", NA)
     leg.col = c("black",c[1],c[2])
     leg.lty = c(NA,1,2)
     leg.pch <- c(19,NA,NA)
     leg.lab[3] <- paste(hdi.params[2]*100, "% Credible Int.", sep = "", collapse = "")

     if(hdi.params[1] == 1){
       hdi.exp <- hdi(x$samples$expectation, credMass = hdi.params[2])
       hdi.alpha <- hdi(x$samples$alpha, credMass = hdi.params[2])
       r.exp<- hdi.exp + 0.5*(hdi.exp-mean.exp)

     }else{
       hdi.alpha <- rep(NA, times = 2)
       leg.lab[3] <-  leg.col[3] <- leg.lty[3] <- NULL
       leg.pch <- leg.pch[1:2]
       hdi.exp <- hdi(x$samples$expectation)
       r.exp <- hdi.exp + 0.5*(hdi.exp-mean.exp)
       hdi.exp <- rep(NA, times = 2)
     }

     d.exp <- density(x$samples$expectation, from = r.exp[1], to = r.exp[2])

     layout(matrix(c(1,1,2,2,3,3,3,3,4,4,4,4), ncol = 4, nrow = 3, byrow = TRUE),heights = c(5,8,2))


     par(mar = c(2,1,2,1))#c(2, 4, 3, 1), cex = 1)
     plot(d.alpha, main = expression(alpha),xlab = "", ylab = "", yaxt = "n",lwd = l.wid,...)
     abline(v = c(-1,1), col = "red", lty = 2,lwd = l.wid)

     par(mar = c(2,1,2,1))#c(2, 1, 3, 2), cex = 1)
     plot(d.exp,xlab = "", ylab = "", xlim = r.exp, main = expression(paste("E[y|",mu,",",alpha,"]")),lwd = l.wid, yaxt = "n",...)
     abline(v = mean.exp, col = c[1],lwd = l.wid)
     abline(v = hdi.exp, col = c[2], lty = 2,lwd = l.wid)

     par(mar = c(4,4,1,1))#c(5, 4, 2, 2),cex = 1)
     plot(0:(length(y)-1),y, type = "p", pch = 19, ylab = "Mass", xlab = "Time Steps", main = "",...)
     abline(h = mean.exp, col = c[1], lwd = l.wid)
     abline(h =hdi.exp, col = c[2], lty = 2, lwd = l.wid)

     par(mar=c(1,1,1,1), xpd = TRUE)#c(1,2,0.5,2), xpd=TRUE)
     plot(1, type = "n", axes=FALSE, xlab="", ylab="", cex = 1,...)
     legend("bottom",horiz = TRUE, legend = leg.lab, col = leg.col,pch = leg.pch, lty = leg.lty, lwd = l.wid)

   }else{
     mean.alpha <- mean(x$samples$alpha)

     r.alpha<- range(x$samples$alpha)
     d.alpha <- density(x$samples$alpha, from = r.alpha[1], to = r.alpha[2])

     layout(matrix(c(1,1,2,2,3,3,3,3), ncol = 4, nrow = 2, byrow = TRUE),heights = c(5,8))


     par(mar = c(2, 4, 3, 1), cex = 1)
     plot(d.alpha, main = expression(alpha),xlab = "", ylab = "", yaxt = "n",lwd = l.wid )
     abline(v = c(-1,1), col = "red", lty = 2,lwd = l.wid)
     if(sum(x$samples$alpha >= 1) > 0){
     x1 <- min(which(d.alpha$x >= 1))
     x2 <- max(which(d.alpha$x <=  r.alpha[2]))
     with(d.alpha, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), density = 50,col="red"))
     }
     if(sum(x$samples$alpha <= -1) > 0){
     x1 <- max(which(d.alpha$x <= -1))
     x2 <- min(which(d.alpha$x >=  r.alpha[1]))
     with(d.alpha, polygon(x=c(x[c(x2,x2:x1,x1)]), y= c(0, y[x2:x1], 0), density = 50,col="red"))
     }

     par(mar = rep(0, times = 4))
     plot(c(0,1),c(0,1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
     samps.out <- signif(mean(x$samples$alpha < 1 & x$samples$alpha > -1)*100, digits = 4)
     samps.out <- paste(samps.out,"%", sep = "")
     text(x = 0.5, y = c(0.6,0.3), c(bquote(.(samps.out) ~ "of the samples of" ~ alpha) , expression("are between (-1,1)")))


     par(mar = c(5, 4, 2, 2),cex = 1)
     plot(0:(length(y)-1),y, type = "p", pch = 19, ylab = "Mass", xlab = "Time Steps", main = "")


   }
  }else if(x$type == "BMB"){

  ## Plots from BMB function
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
    layout(mat = matrix(1:plot.spaces,nrow = nrow.layout, ncol = ncol.layout, byrow = TRUE))
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


}
