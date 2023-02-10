# Gtransition

**Estimate a Stochastic Growth Matrix Based on Length Structure Data**
This package describes a theoretical model expressing the variability
observed in the individual growth, such that each individual in the
population exhibits a growth pattern with a nonlinear trend toward an
expected value. Thus, the growth is represented by the proportion of
individuals in the length class $l$ during a time interval. The
proportion of individuals that grow from length class $l$ to all length
classes $l^{'}$ is represented by a probabilistic density function,
usually gamma distribution or normal distribution; therefore the growth
pattern depends of their parameters, where the mean value indicates the
average growth increment, and the variance explains the individual
variability in growth, consequently both parameters determine the
proportion of individuals going from one length class to another.

## Installation

Install the CRAN version:

``` r
install.packages("Gtransition")
```

Or install de development version:

``` r
# install.packages("devtools")
devtools::install_github("ejosymart/Gtransition")
```

After, that call the package:

``` r
library("Gtransition")
```

## Examples

This is a basic example which shows you how to calculate the transition
growth matrix:

## Mean growth increment (based on von Bertalanffy equation)

``` r
output <- mgi(lowerL = 78, upperL = 202, classL = 4, 
              Linf = 197.42, k = 0.1938, method = "vonB")

output
#> $delta
#>  [1] 20.6867442 19.9820348 19.2773254 18.5726160 17.8679066 17.1631972
#>  [7] 16.4584878 15.7537784 15.0490690 14.3443596 13.6396503 12.9349409
#> [13] 12.2302315 11.5255221 10.8208127 10.1161033  9.4113939  8.7066845
#> [19]  8.0019751  7.2972657  6.5925564  5.8878470  5.1831376  4.4784282
#> [25]  3.7737188  3.0690094  2.3643000  1.6595906  0.9548812  0.2501718
#> [31]  0.0000000
#> 
#> $Laverage
#>  [1]  80  84  88  92  96 100 104 108 112 116 120 124 128 132 136 140 144 148 152
#> [20] 156 160 164 168 172 176 180 184 188 192 196 200
#> 
#> attr(,"class")
#> [1] "Gincrement" "list"
delta    <- output$delta
Laverage <- output$Laverage
```

## Transition growth matrix

``` r
Gmat <- transitionM(lowerL = 78, upperL = 202, classL = 4, 
                   distribution = "gamma", 
                   delta = delta, beta = 0.105, sigma = NULL)
```

# Plots

``` r
plot(Gmat)
```

![](README-unnamed-chunk-4-1.png)<!-- -->

``` r
plot(Gmat, xlab = "XLAB", ylab = "YLAB", adjY = -25,
     col = "grey40", sizeAxis1 = 0.5, sizeAxis2 = 0.5,
     filename = "myplot", 
     savePDF = TRUE, widthPDF = 3, heightPDF = 10, 
     savePNG = TRUE, widthPNG = 300, heightPNG = 1000, resPNG = 110)
```

![](README-unnamed-chunk-5-1.png)<!-- -->

### References

Luquin-Covarrubias M., Morales-Bojorquez E. (2020). Effects of
stochastic growth on population dynamics and management quantities
estimated from an integrated catch-at-length assessment model: *Panopea
globosa* as case study. Ecologial Modeling 440, 109384.
<https://doi.org/10.1016/j.ecolmodel.2020.109384>

Sullivan P.J., Lai H., Galluci V.F. (1990). A Catch-at-Length analysis
that incorporates a stochastic model of growth. Can. J. Fish. Aquat.
Sci. 47: 184-198.
