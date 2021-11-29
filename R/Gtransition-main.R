# Gtransition package: Descripcion corta de lo que hace --------

#' @importFrom stats pnorm pgamma
#' @importFrom utils installed.packages
#' 
#' @title Título de lo que hace.
#'
#' @description Decsripcion de lo que hace
#' @name Gtransition-package
#' @aliases Gtrasition-package Gtrasition
#' @docType package
#' @author Josymar Torrejon-Magallanes <ejosymart@@gmail.com>
#' @details Package: Gtrasition
#' @details Type: Package
#' @details Detalles de lo que hace.
#'
#' In the Gtrasition... 
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
#' @param classL a numeric value ...
#' @param Linf a numeric value ....
#' @param K a numeric value ..
#' @param gm a numeric value .. \code{gm = 1} as default value??.
#' @param dl a numeric value ...  \code{dl = 0.2} as default value??.
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
#' output <- mgi(lowerL = 10, classL = 10, Linf = 60, K = 0.3, method = "vonB")
#'
#' output
#' output$delta
#' output$Laverage
#' @export
mgi <- function(lowerL, classL, Linf,  K, gm = 1, dl = 0.1, growtheq = "vonB"){
  
  if(lowerL >= Linf)
    stop("HEY! 'Linf' must be greather than 'lowerL'")
  
  if(classL > Linf)
    stop("HEY! 'classL' must be lower than 'Linf'")
  
  l_x   <- seq(from = lowerL, to = Linf , by = classL)
  lc_av <- (l_x + classL/2)[-length(l_x)]
  
  estimate <- switch(method,
                     vonB     = .growth_vonB(Linf = Linf, lc_av = lc_av, K = K),
                     Gompertz = .growth_gompertz(Linf = Linf, lc_av = lc_av, K = K),
                     Logistic = .growth_logistic(Linf = Linf, lc_av = lc_av),
                     Schnute  = .growth_schnute(lc_av = lc_av, gm = gm, dl = dl))

  output <- list(delta    = estimate,
                 Laverage = lc_av)
  
  class(output)  <- c("Gincrement", class(output))
  
  return(output)
}




#' Transition Matrix
#'
#' Estimate .....
#' @param lowerL a numeric value that represents...
#' @param classL a numeric value ...
#' @param Linf a numeric value ....
#' @param distribution a character string defining the growth equation to be used.
#' @param delta a numeric vector...
#' @param beta a numeric value...
#' @param sigma a numeric value...
#' @return A list of class 'mtransition'.
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
#' output <- mgi(lowerL = 10, classL = 10, Linf = 60, K = 0.3, method = "vonB")
#' delta <- output$delta
#' 
#' mat <- transitionM(lowerL = 10, classL = 10, Linf = 60, distribution = "gamma", 
#' delta = vonb$delta, beta = 1.5, sigma = NULL)
#' 
#' mat
#' @export
transitionM <- function(lowerL, classL, Linf, distribution = "gamma", 
                        delta, beta = NULL, sigma = NULL){
  
  if(lowerL >= Linf)
    stop("HEY! 'Linf' must be greather than 'lowerL'")
  
  if(classL > Linf)
    stop("HEY! 'classL' must be lower than 'Linf'")
  
  l_x       <- seq(from = lowerL, to = Linf, by = classL)
  lc_av     <- (l_x + classL/2)[-length(l_x)]
  llc       <- length(l_x)
  quantiles <- c(1, l_x+1)
  quantiles <- quantiles[-(length(l_x) + 1)]
  aux       <- matrix(data = 0, nrow = llc, ncol = llc)
  G         <- NULL
  
  # mcdf --------------------------------------------------------------------
  mcdf <- switch(distribution,
                 gamma  = .gamma(llc = llc, quantiles = quantiles, delta = delta, beta = beta),
                 normal = .normal(llc = llc, quantiles = quantiles, delta = delta, sigma = sigma))
  
  rownames(mcdf) <- colnames(mcdf) <- lc_av
  
  # G -----------------------------------------------------------------------
  for(j in seq_len(ncol(mcdf))){
    
    aux <- diff(mcdf[ ,j])/mcdf[nrow(mcdf), j]
    
    G   <- cbind(G,aux)
    
  }
  
  G <- rbind(rep(0, ncol(G)), G)
  rownames(G) <- colnames(G) <- lc_av
  
  out <- list(mcdf = mcdf, G = G)
  
  output <- lapply(out, function(x) round(x, 6))
  
  return(output)
}
