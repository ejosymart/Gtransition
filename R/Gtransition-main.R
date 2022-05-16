# Gtransition package: Estimate a Stochastic Growth Matrix Based on Length Structure data

#' @importFrom stats pnorm pgamma
#' @importFrom utils installed.packages
#' @importFrom graphics axis barplot box mtext par
#' 
#' @title Estimate a stochastic growth matrix based on length structure data
#'
#' @description Describe a theoretical model expressing the variability observed in the individual growth, such that each individual in the population exhibits a growth pattern with a nonlinear trend toward an expected value.
#' @name Gtransition-package
#' @aliases Gtransition-package Gtransition
#' @docType package
#' @author Arelly Ornelas Vargas <aornelasv@@ipn.mx>
#' @author Josymar Torrejon-Magallanes <ejosymart@@gmail.com>
#' @author Marlene Anaid Luquin-Covarrubias <marlene.luquin@@gmail.com>
#' @details Package: Gtrasition
#' @details Type: Package
#' @details The stochastic growth matrix describes a theoretical model expressing the variability observed in the individual growth, such that each individual in the population exhibits a growth pattern with a nonlinear trend toward an expected value.
#' Thus, the growth is represented by the proportion of individuals in the length class \eqn{l} during a time interval. The proportion of individuals that grow from length class \eqn{l} to all length classes \eqn{l^{'}} is represented by a probabilistic density function, usually gamma distribution or normal distribution; therefore, the growth pattern depends on their parameters, where the mean value indicates the average growth increment, and the variance explains the individual variability in growth, consequently both parameters determine the proportion of individuals going from one length class to another.
#'
#' 1) Stochastic growth matrix
#' Individual growth was modeled by using a growth matrix (\eqn{G_{l, l+1}}) which is expressed through a stochastic growth model that defines the probability of each individual growing from one length class to another over a time-step.
#' Mathematically \eqn{G_{l, l+1}} matriz requires estimation of mean growth increments \eqn{\delta_{l}}, which assume length variability from individual to individual estimated by \deqn{\delta_{l} = l_{t+1} - l_{t}}, where \eqn{l_{t+1}} is the length of the individual at time \eqn{t + 1}, and \eqn{l_{t}} is the length of the individual at time \eqn{t}.
#' In this way, The expected mean growth increments were estimated by applying four stochastic growth models.
#' 
#' 1.1) von Bertalanfy stochastic growth model
#' The von Bertalanffy stochastic growth model (VBS) defines an asymptotic curve characterized by an accelerated growth rate in the early stages of development that decreases gradually to attain the asymptotic length (Sullivan et al., 1990; Punt et al., 2010; Cao et al., 2017a; Fisch et al., 2019).
#' 
#' \deqn{\bar{\Delta}_{l} = (L_{\infty} - l_{\ast}) \cdot (1 - e^{-K})}
#' 
#' 
#' 1.2) Gompertz stochastic growth model
#' The re-parameterized Gompertz stochastic growth model (GMS) exhibits an asymmetrical sigmoidal curve with a low inflection point and assumes that growth is not constant throughout the life cycle; thus, younger individuals exhibit faster growth than older individuals (Troynikov et al., 1998; Helidoniotis & Haddon, 2013; Dippold et al., 2017).
#' 
#' \deqn{\bar{\Delta}_{l} = L_{\infty} \cdot \frac{l_{\ast}}{L_{\infty}}^{exp(-K)} - l_{\ast}}
#' 
#' 
#' 1.3) Logistic stochastic growth model
#' The Logistic stochastic growth model (LGS) describes a symmetrical sigmoidal curve and denotes several growth possibilities, which can be spread to maximum lengths, allowing the description of both determinate and indeterminate growth (Haddon et al., 2008; Helidoniotis et al., 2011).
#' 
#' \deqn{\bar{\Delta}_{l} = \frac{Max\Delta_{l}}{1 + exp(-LN(19) \cdot (\frac{(l_{\ast})-L_{50}}{L_{95} - L{50}}))}}
#' 
#' 
#' 1.4) Schnute stochastic growth model
#' The Schnute stochastic growth matrix (SCS) was used assuming the parameters \eqn{\delta} \eqn{\neq} 0, \eqn{\gamma} \eqn{\neq} 0,
#' where \eqn{\delta} represents a constant relative rate of relative growth rate, and \eqn{\gamma} is the incremental relative rate of relative growth rate (Schnute, 1981).
#' SCS SCS is a general growth model with high flexibility describing a variety of growth patterns (asymptotic, linear, and exponential), including properties such as growth acceleration, 
#' asymptotic limits and inflection points, depending on the parameter values (Baker et al., 1991). According to Schnute (1981), if \eqn{\gamma} = 1,
#' then the model describes an asymptotic growth pattern that corresponds to the von Bertalanffy shape>
#' 
#' \deqn{\bar{\Delta}_{l} = -l_{\ast} + (l^{\gamma}_{\ast} \cdot exp^{-\gamma} + L^{\gamma}_{\infty} \cdot (1 - exp^{-\gamma}))^{\frac{1}{\gamma}}}
#' 
#' where \eqn{\bar{\Delta}_{l}} is the expected mean growth increment for length class \eqn{l}, \eqn{l_{\ast}} represents the midlength of the length class \eqn{l},
#' \eqn{L_{\infty}} is the asymptotic length where the mean growth increment is zero (VBS, GMS and SCS), \eqn{K} represents the growth rate (VBS and GMS), 
#' Max\Delta_{l} is the maximum growth increment, \eqn{L_{50}} is the initial length that produces a growth increment of 0.5 times Max\Delta_{l}, and
#' \eqn{L_{95}} is the initial length at 0.05 times Max\Delta_{l} (LGS). Te equation LGS uses \eqn{-LN(19)}, thus expressing a Logistic curve; if LN(19) is used then an inverse 
#' logistic curve could be modeled (Baker et al., 1991; Haddon et al., 2008; Helidoniotis et al., 2011).
#' 
#' 
#' @references Luquin-Covarrubias M., Morales-Bojorquez E. (2020). Effects of stochastic growth on population dynamics and management quantities estimated from an integrated catch-at-length assessment model: Panopea globosa as case study. Ecologial Modeling 440, 109384. https://doi.org/10.1016/j.ecolmodel.2020.109384
#' @references Sullivan P.J., Lai H., Galluci V.F. (1990). A Catch-at-Length analysis that incorporates a stochastic model of growth. Can. J. Fish. Aquat. Sci. 47: 184-198.
#' @concept Gtransition
#' @concept growth
#' @concept transition
#' @concept distributions
#' @examples
#' #See examples for functions.

NULL
#' Mean growth increment
#'
#' Estimate mean growth increment, for the individuals in length class I is then the average change in length of individuals initially in length class
#' @param lowerL a numeric value that represents...
#' @param upperL a numeric value that represents...
#' @param classL a numeric value ...
#' @param Linf a numeric value ....
#' @param K a numeric value ..
#' @param gm a numeric value .. \code{gm = 1} as default value??.
#' @param dl a numeric value ...  \code{dl = 0.1} as default value??.
#' @param method a character string defining the growth equation to be used.
#' @return A list of class 'Gincrement'.
#'
#' \code{delta} the mean growth increment
#'
#' \code{Laverage} the average length.
#'
#' @details Estimate mean growth increment
#'
#'
#' @references Sullivan P.J., Lai H., Galluci V.F. (1990). A Catch-at-Length analysis that incorporates a stochastic model of growth. Can. J. Fish. Aquat. Sci. 47: 184-198.
#' @examples
#' output <- mgi(lowerL = 78, upperL = 202, classL = 4, Linf = 197.42, K = 0.1938, method = "vonB")
#'
#' output
#' output$delta
#' output$Laverage
#' @export
mgi <- function(lowerL, upperL, classL, Linf,  K, gm = 1, dl = 0.1, method = "vonB"){
  
  if(lowerL >= Linf)
    stop("HEY! 'Linf' must be greather than 'lowerL'")
  
  if(lowerL >= upperL)
    stop("HEY! 'upperL' must be greather than 'lowerL'")
  
  if(classL > Linf)
    stop("HEY! 'classL' must be lower than 'Linf'")
  
  l_x   <- seq(from = lowerL, to = upperL, by = classL)
  lc_av <- (l_x + classL/2)[-length(l_x)]
  
  estimate <- switch(method,
                     vonB     = .growth_vonB(Linf = Linf, lc_av = lc_av, K = K),
                     Gompertz = .growth_gompertz(Linf = Linf, lc_av = lc_av, K = K),
                     Logistic = .growth_logistic(Linf = Linf, lc_av = lc_av),
                     Schnute  = .growth_schnute(lc_av = lc_av, gm = gm, dl = dl, Linf = Linf))
  
  output <- list(delta    = estimate,
                 Laverage = lc_av)
  
  class(output)  <- c("Gincrement", class(output))
  
  return(output)
}




#' Transition Matrix
#'
#' Estimate .....
#' @param lowerL a numeric value that represents...
#' @param upperL a numeric value that represents...
#' @param classL a numeric value ...
#' @param distribution a character string defining the growth equation to be used.
#' @param delta a numeric vector...
#' @param beta a numeric value...
#' @param sigma a numeric value...
#' @return A list of class 'Mtransition'.
#'
#' \code{mcdf} the mean growth increment
#'
#' \code{G} the average length.
#'
#' @details Estimate .....
#'
#'
#' @references Sullivan P.J., Lai H., Galluci V.F. (1990). A Catch-at-Length analysis that incorporates a stochastic model of growth. Can. J. Fish. Aquat. Sci. 47: 184-198.
#' @examples
#' output <- mgi(lowerL = 78, upperL = 202, classL = 4, Linf = 197.42, K = 0.1938, method = "vonB")
#' delta <- output$delta
#' 
#' mat <- transitionM(lowerL = 78, upperL = 202, classL = 4, distribution = "gamma", 
#' delta = delta, beta = 0.105, sigma = NULL)
#' 
#' mat
#' @export
transitionM <- function(lowerL, upperL, classL, distribution = "gamma", 
                        delta, beta = NULL, sigma = NULL){
  
  if(lowerL >= upperL)
    stop("HEY! 'upperL' must be greather than 'lowerL'")
  
  l_x       <- seq(from = lowerL, to = upperL, by = classL)
  lc_av     <- (l_x + classL/2)[-length(l_x)]
  llc       <- length(lc_av)
  aux       <- matrix(data = 0L, nrow = llc, ncol = llc)
  G         <- NULL
  
  # mcdf --------------------------------------------------------------------
  mcdf <- switch(distribution,
                 gamma  = .gamma(llc = llc, lc_av = lc_av, delta = delta, beta = beta),
                 normal = .normal(llc = llc, lc_av = lc_av, delta = delta, sigma = sigma))
  
  rownames(mcdf) <- colnames(mcdf) <- lc_av
  
  # G -----------------------------------------------------------------------
  for(j in seq_len(ncol(mcdf))){
    
    aux <- diff(mcdf[ ,j])/mcdf[nrow(mcdf), j]
    
    G   <- cbind(G,aux)
    
  }
  
  G <- rbind(rep(0, ncol(G)), G)
  rownames(G) <- colnames(G) <- lc_av
  
  output <- G
  
  class(output)  <- c("Mtransition", class(output))
  
  return(output)
}



#' Plot method for Mtransition class
#'
#' @param x object of class 'Mtransition'.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param sizeAxis1 a number for the axis X.
#' @param sizeAxis2 a number for the axis Y.
#' @param adjY a number for y label position 
#' @param col color for the barplot.
#' @param \dots Additional arguments to the plot method.
#' @examples
#' output <- mgi(lowerL = 78, upperL = 202, classL = 4, Linf = 197.42, K = 0.1938, method = "vonB")
#' delta <- output$delta
#' 
#' mat <- transitionM(lowerL = 78, upperL = 202, classL = 4, distribution = "gamma", 
#' delta = delta, beta = 0.105, sigma = NULL)
#' 
#' plot(mat)
#' @export
#' @method plot Mtransition
plot.Mtransition <- function(x, xlab = "X-Text", ylab = "Y-Text", col = "grey45", 
                             sizeAxis1 = 0.85, sizeAxis2 = 0.5, adjY = -15.5, ...){
  
  if (!inherits(x, "Mtransition"))
    stop("Use only with 'Mtransition' objects")
  
  data     <- x
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(ncol(data), 1), 
      mai   = c(0.05, 0.4, 0.05, 0.2), 
      oma   = c(4, 3, 0, 4)) 
  
  for(i in rev(seq_len(ncol(data)))){
    barplot(data[,i], names.arg = rownames(data), 
            ylim = c(0, 1.1*max(data)), 
            las = 1, 
            space = 0, 
            border = NA, 
            col = col, 
            xaxt = "n", 
            yaxt = "n")
    box()
    axis(side = 2, at = c(0, 0.5), las = 2, cex.axis = sizeAxis2)
  }
  axis(side = 1, at = seq_along(rownames(data)), 
       labels = rownames(data), las = 2, 
       cex.axis = sizeAxis1)
  mtext(text = xlab, side = 1, line = 2.75)
  mtext(text = ylab, side = 2, line = 2.75, adj = adjY)
  
  
  return(invisible(NULL))
}
