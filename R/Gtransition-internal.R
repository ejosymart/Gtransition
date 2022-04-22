.growth_vonB <- function(Linf, lc_av, K){
  
  d <- (Linf - lc_av)*(1- exp(-K))
  
  return(d)
}

.growth_gompertz <- function(Linf, lc_av, K){
  
  d <- (Linf*(lc_av/Linf)^exp(-K))-lc_av
  
  return(d)
}

.growth_logistic <- function(Linf, lc_av){
  
  dmx <- Linf - lc_av[1]
  
  d <- dmx/(1+exp(-log(19)*(lc_av - dmx*0.5)/(dmx*0.05 - dmx*0.5)))
  
  return(d)
}


.growth_schnute <- function(lc_av, gm, dl, Linf){
  
  d <- -lc_av + (lc_av^gm * exp(-dl) + Linf^gm*(1 - exp(-dl)))^(1/gm)
  
  return(d)
}



.gamma <- function(llc, lc_av, delta, beta){
  
  aux  <- matrix(data = 0L, nrow = llc, ncol = llc)
  
  mcdf <- NULL
  for(j in seq_len(llc)){
    alpha <- (delta + lc_av)/beta
    aux[, j] <- pgamma(q = lc_av, shape = alpha[j], scale = beta, lower.tail = T)

    mcdf <- cbind(mcdf, aux)
  }
  
  mcdf[upper.tri(mcdf)] <- 0
  
  rownames(mcdf) <- colnames(mcdf) <-lc_av
  
  return(mcdf)
}


.normal <- function(llc, lc_av, delta, sigma){
  
  aux  <- matrix(data = 0L, nrow = llc, ncol = llc)
  
  mcdf <- NULL
  for(j in seq_len(llc)){
    aux[, j] <- pnorm(q = lc_av, mean = delta[j], sd = sigma, lower.tail = T)
    
    mcdf <- cbind(mcdf, aux)
  }
  
  mcdf[upper.tri(mcdf)] <- 0
  
  rownames(mcdf) <- colnames(mcdf) <-lc_av
  
  return(mcdf)
}