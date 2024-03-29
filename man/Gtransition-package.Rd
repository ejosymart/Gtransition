% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gtransition-main.R
\docType{package}
\name{Gtransition-package}
\alias{Gtransition-package}
\alias{Gtransition}
\title{Estimate a stochastic growth matrix based on length structure data}
\description{
Describe a theoretical model expressing the variability observed in the individual growth, such that each individual in the population exhibits a growth pattern with a nonlinear trend toward an expected value. It includes estimation of mean growth increment \eqn{\bar \Delta_l} for the individuals in length class \eqn{l} and the probability of each individual growing from one length class to another over a time-step \eqn{G_{l, l+1}} based on two probabilistic density functions: gamma and normal distributions, and some basic plots.
}
\details{
Package: Gtransition

Type: Package

The stochastic growth matrix describes a theoretical model expressing the variability observed in the individual growth, 
such that each individual in the population exhibits a growth pattern with a nonlinear trend toward an expected value.
Thus, the growth is represented by the proportion of individuals in the length class \eqn{l} during a time interval. 
The proportion of individuals that grow from length class \eqn{l} to all length classes \eqn{l^{'}} is represented 
by a probabilistic density function, usually gamma distribution or normal distribution; therefore, the growth pattern 
depends on their parameters, where the mean value indicates the average growth increment, and the variance explains 
the individual variability in growth, consequently both parameters determine the proportion of individuals going 
from one length class to another.

\describe{

\item{Stochastic growth matrix:}{}

Individual growth was modeled by using a growth matrix (\eqn{G_{l, l+1}}) which is expressed through a stochastic growth model 
that defines the probability of each individual growing from one length class to another over a time-step.
Mathematically \eqn{G_{l, l+1}} matrix requires estimation of mean growth increments \eqn{\delta_{l}}, which assume length 
variability from individual to individual estimated by \eqn{\delta_{l} = l_{t+1} - l_{t}}, where \eqn{l_{t+1}} is the length 
of the individual at time \eqn{t + 1}, and \eqn{l_{t}} is the length of the individual at time \eqn{t}.
In this way, the expected mean growth increments were estimated by applying four stochastic growth models.

\enumerate{
\item von Bertalanfy stochastic growth model (VBS):

The VBS growth model defines an asymptotic curve characterized by an accelerated growth 
rate in the early stages of development that decreases gradually to attain the asymptotic length 
(Sullivan et al., 1990; Punt et al., 2010; Cao et al., 2017a; Fisch et al., 2019).

\deqn{\bar{\Delta}_{l} = (L_{\infty} - l_{\ast}) (1 - e^{-k})}


\item Gompertz stochastic growth model (GMS):

The re-parameterized Gompertz stochastic growth model (GMS) exhibits an asymmetrical sigmoidal curve with a low inflection 
point and assumes that growth is not constant throughout the life cycle; thus, younger individuals exhibit faster growth 
than older individuals (Troynikov et al., 1998; Helidoniotis & Haddon, 2013; Dippold et al., 2017).

\deqn{\bar{\Delta}_{l} = L_{\infty} \left(\frac{l_{\ast}}{L_{\infty}}\right)^{\exp(-k)} - l_{\ast}}


\item Logistic stochastic growth model:

The logistic stochastic growth model (LGS) describes a symmetrical sigmoidal curve and denotes several growth possibilities, 
which can be spread to maximum lengths, allowing the description of both determinate and indeterminate growth 
(Haddon et al., 2008; Helidoniotis et al., 2011).

\deqn{\bar{\Delta}_{l} = \frac{\textup{Max} \; \Delta_{l}}{1 + \exp \left(-\ln(19) (\frac{l_{\ast}-L_{50}}{L_{95} - L{50}})\right)}}


\item Schnute stochastic growth model:

The Schnute stochastic growth matrix (SCS) was used assuming the parameters \eqn{\delta} \eqn{\neq} 0, \eqn{\gamma} \eqn{\neq} 0,
where \eqn{\delta} represents a constant relative rate of relative growth rate, and \eqn{\gamma} is the incremental relative rate 
of relative growth rate (Schnute, 1981).
SCS is a general growth model with high flexibility describing a variety of growth patterns (asymptotic, linear, and exponential), 
including properties such as growth acceleration, asymptotic limits and inflection points, depending on the parameter values 
(Baker et al., 1991). According to Schnute (1981), if \eqn{\gamma} = 1, then the model describes an asymptotic growth pattern 
that corresponds to the von Bertalanffy shape:

\deqn{\bar{\Delta}_{l} = -l_{\ast} + (l^{\gamma}_{\ast} \exp^{-\delta} + L^{\gamma}_{\infty} (1 - \exp^{-\delta}))^{\frac{1}{\gamma}},}

where \eqn{\bar{\Delta}_{l}} is the expected mean growth increment for length class \eqn{l}, 
\eqn{l_{\ast}} represents the midlength of the length class \eqn{l},
\eqn{L_{\infty}} is the asymptotic length where the mean growth increment is zero (VBS, GMS and SCS), 
\eqn{k} represents the growth rate (VBS and GMS), 
\eqn{\textup{Max}\;\Delta_{l}} is the maximum growth increment, 
\eqn{L_{50}} is the initial length that produces a growth increment of \eqn{0.5} times \eqn{\textup{Max}\;\Delta_{l}}, 
and \eqn{L_{95}} is the initial length at \eqn{0.05} times \eqn{\textup{Max}\;\Delta_{l}} (LGS). 
The equation LGS uses \eqn{-\ln(19)}, thus expressing a logistic curve; if if \eqn{\ln(19)} is used then an inverse 
logistic curve could be modeled (Baker et al., 1991; Haddon et al., 2008; Helidoniotis et al., 2011).
}

\item{Variability:}{}
To describe the variability in the mean growth increments, the \eqn{G_{l, l+1}} matrix was based 
on two probabilistic density functions: gamma and normal distributions. 
Both functions define the probability region where individuals may grow, including the probability 
that the increment in length does not occur and the individuals remain in their original 
length class (Haddon, 2011). Thus, the probabilities of growth increments were estimated:

\enumerate{
\item Assuming a gamma distribution:

\deqn{g(\Delta_{l}|\alpha_{l}\beta_{g}) = \frac{1}{\beta^{\alpha_{l}}_{g} \Gamma(\alpha_{l})} \Delta_{l}^{\alpha_{l}-1} e^{-\frac{\Delta_{l}}{\beta_{g}}},}

where \eqn{\alpha_{l}} is the scale parameter, \eqn{\beta_{g}} is the shape parameter 
and \eqn{\Gamma} is the gamma funtion for the \eqn{\alpha_{l}} parameter.
The mean change in length is \eqn{\bar{\Delta}_{l} = \alpha_{l} \beta_{g}} 
and the variance is \eqn{\sigma_{\Gamma}^{2} = \alpha_{l} \beta_{g}^{2}}.
The expected proportion of individuals growing from length class \eqn{l} 
to length class \eqn{l + 1} can be found by integrating over the length range 
\eqn{l + 1_{1}, l + 1_{2}} which represent the lower and upper ends of length classes, 
respectively (Quinn & Deriso, 1999; Haddon, 2011):

\deqn{G_{l, l+1} = \int_{l+1_{1}}^{l+1_{2}} g(x|\alpha_{1}, \beta_{g}) dx}



\item Assuming a normal distribution:

\deqn{X_{k}= \int_{L_{j}}^{L_{j+1}} \frac{1}{\sigma_{k}\sqrt{2\pi}} \exp\left(-\frac{L-(\tilde{L}_{i} + I_{k})}{2(\sigma^2_k}\right) dL,}

where \eqn{\sigma_{k}} determines the variability in growth increment (k) for individuals, \eqn{\tilde{L}_{i}} is the midpoint of the length class \eqn{i},
\eqn{I_{k}} is the growth increment. The normal distribution defines the probability that an individual in length class \eqn{i} grows into size class \eqn{j}
during each time-step (Haddon 2011).
}
}
}
\examples{
#See examples for functions mgi() and transitionM.
}
\references{
Baker, T.T., Lafferty, R., Quinn II, T.J., 1991. A general growth model for mark-recapture data. Fisheries Research. 11(3-4), 257-281. https://doi.org/10.1016/0165-7836(91)90005-Z.

Cao, J., Chen, Y., Richards, R.A., 2017a. Improving assessment of Pandalus stocks using a seasonal, size-structured assessment model with environmental variables. Part I: Model description and application. Canadian Journal of Fisheries and Aquatic Sciences, 74(3), 349-362. https://doi.org/10.1139/cjfas-2016-0020.

Dippold, D.A., Leaf, R.T., Franks, J.S., Hendon, J.R., 2017. Growth, mortality, and movement of cobia (\eqn{Rachycentron} \eqn{canadum}). Fishery Bulletin, 115(4). doi: 10.7755/FB.115.4.3

Fisch, N.C., Bence, J.R., Myers, J.T., Berglund, E.K., Yule, D.L., 2019. A comparison of age-and size-structured assessment models applied to a stock of cisco in Thunder Bay, Ontario. Fisheries Research, 209, 86-100. https://doi.org/10.1016/j.fishres.2018.09.014.

Haddon, M., Mundy, C., Tarbath, D., 2008. Using an inverse-logistic model to describe growth increments of blacklip abalone (Haliotis rubra) in Tasmania. Fishery Bulletin, 106(1), 58-71.

Haddon, M., 2011. Modelling and quantitative methods in fisheries. CRC press. Second ed. Boca Raton, London, New York.

Helidoniotis, F., Haddon, M., Tuck, G., Tarbath, D., 2011. The relative suitability of the von Bertalanffy, Gompertz and inverse logistic models for describing growth in blacklip abalone populations (Haliotis rubra) in Tasmania, Australia. Fisheries Research, 112(1-2), 13-21. https://doi.org/10.1016/j.fishres.2011.08.005.

Helidoniotis, F., Haddon, M., 2013. Growth models for fisheries: The effect of unbalanced sampling error on model selection, parameter estimation, and biological predictions. Journal of Shellfish Research, 32(1), 223-236. https://doi.org/10.2983/035.032.0129.

Schnute, J., 1981. A versatile growth model with statistically stable parameters. Canadian Journal of Fisheries and Aquatic Sciences, 38(9), 1128-1140. https://doi.org/10.1139/f81-153.

Sullivan, P.J., Lai, H.L., Gallucci, V.F., 1990. A catch-at-length analysis that incorporates a stochastic model of growth. Canadian Journal of Fisheries and Aquatic Sciences, 47(1), 184-198. https://doi.org/10.1139/f90-021.

Punt, A.E., Deng, R.A., Dichmont, C.M., Kompas, T., Venables,W.N., Zhou, S., Pascoe, S., Hutton, T., Kenyon, R., van der Velde, T., Kienzle, M., 2010. Integrating size-structured assessment and bioeconomic management advice in Australia's northern prawn fishery. ICES Journal of Marine Science, 67(8), 1785-1801. https://doi.org/10.1093/icesjms/fsq037.

Quinn II, T.J., Deriso, R.B., 1999. Quantitative fish dynamics. First ed. New York, Oxford.

Troynikov, V.S., Day, R.W., Leorke, A.M., 1998. Estimation of seasonal growth parameters using a stochastic Gompertz model for tagging data. Journal of Shellfish Research, 17, 833-838.
}
\author{
Arelly Ornelas-Vargas <aornelasv@ipn.mx>

Josymar Torrejon-Magallanes <ejosymart@gmail.com>

Marlene Anaid Luquin-Covarrubias <marlene.luquin@gmail.com>
}
\concept{Gtransition}
\concept{distribution}
\concept{growth}
\concept{transition}
