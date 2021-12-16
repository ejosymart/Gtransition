# Gtransition

**Estimate growth trastion matrix** This package allows to…

## Installation

Install the CRAN version:

Or install de development version:

After, that call the package:

## Examples

This is a basic example which shows you how to calculate the transition
growth matrix:

## Mean growth increment (based on von Bertalanffy equation)

``` r
output <- mgi(lowerL = 10, classL = 10, Linf = 60, K = 0.3, method = "vonB")

output
#> $delta
#> [1] 11.663180  9.071362  6.479544  3.887727  1.295909
#> 
#> $Laverage
#> [1] 15 25 35 45 55
#> 
#> attr(,"class")
#> [1] "Gincrement" "list"
output$delta
#> [1] 11.663180  9.071362  6.479544  3.887727  1.295909
output$Laverage
#> [1] 15 25 35 45 55
```

## Transition growth matrix

``` r
mat <- transitionM(lowerL = 10, classL = 10, Linf = 60, 
                   distribution = "gamma", delta = output$delta, 
                   beta = 1.5, sigma = NULL)

mat
#> $mcdf
#>          15       25       35       45       55
#> 15 0.000001 0.000000 0.000000 0.000000 0.000000
#> 25 0.483565 0.000062 0.000000 0.000000 0.000000
#> 35 0.973265 0.733562 0.002585 0.000000 0.000000
#> 45 0.999612 0.994192 0.913049 0.059277 0.000000
#> 55 0.999997 0.999954 0.999253 0.986472 0.554959
#> 
#> $G
#>          15       25       35      45 55
#> 15 0.000000 0.000000 0.000000 0.00000  0
#> 25 0.483566 0.000062 0.000000 0.00000  0
#> 35 0.489701 0.733534 0.002587 0.00000  0
#> 45 0.026347 0.260642 0.911144 0.06009  0
#> 55 0.000385 0.005762 0.086269 0.93991  1
```

### References

Luquin-Covarrubias M., Morales-Bojorquez E. (2020). Effects of
stochastic growth on population dynamics and management quantities
estimated from an integrated catch-at-length assessment model: Panopea
globosa as case study. Ecologial Modeling 440, 109384.
<https://doi.org/10.1016/j.ecolmodel.2020.109384>

Sullivan P.J., Lai H., Galluci V.F. (1990). A Catch-at-Length analysis
that incorporates a stochastic model of growth. Can. J. Fish. Aquat.
Sci. 47: 184-198.
