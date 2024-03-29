% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gtransition-main.R
\name{transitionM}
\alias{transitionM}
\title{Transition Matrix}
\usage{
transitionM(
  lowerL,
  upperL,
  classL,
  distribution = "gamma",
  delta,
  beta = NULL,
  sigma = NULL
)
}
\arguments{
\item{lowerL}{a numeric value that represents the smallest observed  size.}

\item{upperL}{a numeric value that represents the highest observed size.}

\item{classL}{a numeric value that represents the range length classes.}

\item{distribution}{a character string defining the growth equation to be used. One of "gamma" or "normal".}

\item{delta}{a numeric vector that represents mean growth increments \eqn{\bar \Delta_l}. The output of function mgi().}

\item{beta}{a numeric value that represents  the shape parameter of gamma distribution density function. Required when "gamma" distribution is selected.}

\item{sigma}{a numeric value the represents the variability  of normal distribution density function. Required when "normal" distribution is selected.}
}
\value{
A list of class 'Mtransition'.

\code{mcdf} a matrix that contain the probabilities  of growth increments.

\code{G} a matrix that contain the expected proportion of individuals growing from class \eqn{l} to  length class \eqn{l +1}.
}
\description{
Estimates the probability of each individual growing from one length class to another over a time-step \eqn{G_{l, l+1}} based on two probabilistic density functions: gamma and normal distributions.
}
\details{
The probabilistic density function defines the probability region where individuals 
may grow considering the probability that the increment in length does not occur and the individuals 
remain in their original length class (Haddon, 2011).
}
\examples{
output <- mgi(lowerL = 78, upperL = 202, classL = 4, Linf = 197.42, k = 0.1938, method = "vonB")
delta <- output$delta

mat <- transitionM(lowerL = 78, upperL = 202, classL = 4, distribution = "gamma", 
delta = delta, beta = 0.105, sigma = NULL)

mat
}
\references{
Luquin-Covarrubias, M.A. and E. Morales-Bojórquez. 2021. Effects of stochastic growth on population dynamics and management quantities estimated from an integrated catch-at-length assessment model: Panopea globosa as case study. Ecological Modelling. 440: 109384. https://doi.org/10.1016/j.ecolmodel.2020.109384.
}
