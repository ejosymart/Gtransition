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
output <- mgi(lowerL = 78, upperL = 202, classL = 4, 
              Linf = 197.42, K = 0.1938, method = "vonB")

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
mat <- transitionM(lowerL = 78, upperL = 202, classL = 4, 
                   distribution = "gamma", 
                   delta = delta, beta = 0.105, sigma = NULL)
 
mat$mcdf
#>           80       84       88       92       96      100      104      108
#> 80  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 84  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 88  0.000024 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 92  0.003056 0.000083 0.000001 0.000000 0.000000 0.000000 0.000000 0.000000
#> 96  0.073028 0.006760 0.000255 0.000004 0.000000 0.000000 0.000000 0.000000
#> 100 0.420380 0.113095 0.013591 0.000694 0.000015 0.000000 0.000000 0.000000
#> 104 0.845843 0.506394 0.164491 0.025102 0.001692 0.000049 0.000001 0.000000
#> 108 0.986366 0.887098 0.589114 0.226428 0.042997 0.003749 0.000146 0.000003
#> 112 0.999621 0.991282 0.918854 0.665450 0.296999 0.068871 0.007622 0.000389
#> 116 0.999997 0.999781 0.994490 0.942686 0.733328 0.373441 0.103910 0.014339
#> 120 1.000000 0.999998 0.999874 0.996554 0.960168 0.791695 0.452544 0.148629
#> 124 1.000000 1.000000 0.999999 0.999928 0.997866 0.972728 0.840382 0.531073
#> 128 1.000000 1.000000 1.000000 0.999999 0.999959 0.998690 0.981584 0.879888
#> 132 1.000000 1.000000 1.000000 1.000000 1.000000 0.999977 0.999202 0.987723
#> 136 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 0.999987 0.999518
#> 140 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 0.999993
#> 144 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 148 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 152 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 156 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 160 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 164 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 168 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 172 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 176 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 180 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 184 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 188 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 192 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 196 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 200 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#>          112      116      120      124      128      132      136      140
#> 80  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 84  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 88  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 92  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 96  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 100 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 104 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 108 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 112 0.000009 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 116 0.000942 0.000029 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 120 0.025160 0.002090 0.000084 0.000000 0.000000 0.000000 0.000000 0.000000
#> 124 0.202710 0.041452 0.004287 0.000219 0.000000 0.000000 0.000000 0.000000
#> 128 0.606145 0.264981 0.064524 0.008189 0.000527 0.000000 0.000000 0.000000
#> 132 0.911148 0.675483 0.333539 0.095408 0.014661 0.001170 0.000000 0.000000
#> 136 0.991912 0.935319 0.737535 0.405982 0.134674 0.024739 0.002418 0.000000
#> 140 0.999711 0.994730 0.953620 0.791486 0.479699 0.182290 0.039553 0.004678
#> 144 0.999996 0.999828 0.996601 0.967211 0.837161 0.552153 0.237565 0.060199
#> 148 1.000000 0.999998 0.999898 0.997829 0.977125 0.874895 0.621124 0.299201
#> 152 1.000000 1.000000 0.999999 0.999940 0.998625 0.984239 0.905371 0.684866
#> 156 1.000000 1.000000 1.000000 0.999999 0.999965 0.999137 0.989267 0.929474
#> 160 1.000000 1.000000 1.000000 1.000000 1.000000 0.999979 0.999462 0.992771
#> 164 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 0.999988 0.999667
#> 168 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 0.999993
#> 172 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 176 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 180 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 184 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 188 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 192 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 196 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#> 200 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#>          144      148      152      156      160      164      168      172
#> 80  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 84  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 88  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 92  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 96  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 100 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 104 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 108 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 112 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 116 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 120 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 124 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 128 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 132 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 136 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 140 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 144 0.008518 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 148 0.087591 0.014672 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 152 0.365416 0.122315 0.024008 0.000000 0.000000 0.000000 0.000000 0.000000
#> 156 0.742183 0.434139 0.164512 0.037476 0.000000 0.000000 0.000000 0.000000
#> 160 0.948171 0.792430 0.503214 0.213818 0.056012 0.000000 0.000000 0.000000
#> 164 0.995181 0.962415 0.835449 0.570606 0.269362 0.080430 0.000000 0.000000
#> 168 0.999796 0.996819 0.973086 0.871480 0.634553 0.329840 0.111311 0.000000
#> 172 0.999996 0.999875 0.997919 0.980956 0.901045 0.693681 0.393629 0.148903
#> 176 1.000000 0.999998 0.999924 0.998650 0.986676 0.924844 0.747044 0.458938
#> 180 1.000000 1.000000 0.999999 0.999954 0.999131 0.990778 0.943660 0.794125
#> 184 1.000000 1.000000 1.000000 0.999999 0.999973 0.999445 0.993681 0.958290
#> 188 1.000000 1.000000 1.000000 1.000000 1.000000 0.999984 0.999648 0.995712
#> 192 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 0.999990 0.999779
#> 196 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 0.999994
#> 200 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
#>          176      180      184      188      192      196      200
#> 80  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 84  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 88  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 92  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 96  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 100 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 104 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 108 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 112 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 116 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 120 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 124 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 128 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 132 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 136 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 140 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 144 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 148 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 152 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 156 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 160 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 164 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 168 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 172 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 176 0.193062 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 180 0.523968 0.243229 0.000000 0.000000 0.000000 0.000000 0.000000
#> 184 0.834790 0.587051 0.298462 0.000000 0.000000 0.000000 0.000000
#> 188 0.969487 0.869220 0.646763 0.357505 0.000000 0.000000 0.000000
#> 192 0.997116 0.977930 0.897829 0.701996 0.418898 0.000000 0.000000
#> 196 0.999862 0.998077 0.984210 0.921187 0.751989 0.481088 0.000000
#> 200 0.999997 0.999914 0.998729 0.988819 0.939942 0.796322 0.503047

mat$G
#>           80       84       88       92       96      100      104      108
#> 80  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 84  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 88  0.000024 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 92  0.003032 0.000083 0.000001 0.000000 0.000000 0.000000 0.000000 0.000000
#> 96  0.069972 0.006677 0.000254 0.000004 0.000000 0.000000 0.000000 0.000000
#> 100 0.347352 0.106335 0.013336 0.000690 0.000015 0.000000 0.000000 0.000000
#> 104 0.425463 0.393300 0.150900 0.024408 0.001677 0.000049 0.000001 0.000000
#> 108 0.140522 0.380703 0.424623 0.201327 0.041305 0.003700 0.000146 0.000003
#> 112 0.013255 0.104184 0.329740 0.439022 0.254001 0.065122 0.007476 0.000387
#> 116 0.000376 0.008499 0.075636 0.277236 0.436330 0.304570 0.096288 0.013950
#> 120 0.000003 0.000217 0.005385 0.053868 0.226839 0.418254 0.348633 0.134290
#> 124 0.000000 0.000002 0.000125 0.003374 0.037698 0.181032 0.387839 0.382444
#> 128 0.000000 0.000000 0.000001 0.000071 0.002093 0.025962 0.141202 0.348815
#> 132 0.000000 0.000000 0.000000 0.000001 0.000041 0.001287 0.017619 0.107835
#> 136 0.000000 0.000000 0.000000 0.000000 0.000000 0.000023 0.000784 0.011795
#> 140 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000013 0.000475
#> 144 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000007
#> 148 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 152 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 156 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 160 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 164 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 168 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 172 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 176 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 180 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 184 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 188 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 192 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 196 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 200 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#>          112      116      120      124      128      132      136      140
#> 80  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 84  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 88  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 92  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 96  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 100 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 104 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 108 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 112 0.000009 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 116 0.000933 0.000029 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 120 0.024218 0.002061 0.000084 0.000000 0.000000 0.000000 0.000000 0.000000
#> 124 0.177551 0.039362 0.004204 0.000219 0.000000 0.000000 0.000000 0.000000
#> 128 0.403435 0.223529 0.060236 0.007970 0.000527 0.000000 0.000000 0.000000
#> 132 0.305002 0.410502 0.269015 0.087218 0.014134 0.001170 0.000000 0.000000
#> 136 0.080764 0.259837 0.403996 0.310574 0.120014 0.023569 0.002418 0.000000
#> 140 0.007799 0.059411 0.216085 0.385504 0.345024 0.157551 0.037135 0.004678
#> 144 0.000285 0.005097 0.042981 0.175725 0.357463 0.369863 0.198012 0.055521
#> 148 0.000004 0.000170 0.003296 0.030618 0.139964 0.322742 0.383559 0.239002
#> 152 0.000000 0.000002 0.000101 0.002111 0.021500 0.109344 0.284246 0.385665
#> 156 0.000000 0.000000 0.000001 0.000060 0.001339 0.014898 0.083897 0.244608
#> 160 0.000000 0.000000 0.000000 0.000001 0.000035 0.000842 0.010195 0.063297
#> 164 0.000000 0.000000 0.000000 0.000000 0.000000 0.000020 0.000526 0.006896
#> 168 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000012 0.000326
#> 172 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000007
#> 176 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 180 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 184 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 188 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 192 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 196 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 200 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#>          144      148      152      156      160      164      168      172
#> 80  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 84  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 88  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 92  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 96  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 100 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 104 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 108 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 112 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 116 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 120 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 124 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 128 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 132 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 136 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 140 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 144 0.008518 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 148 0.079073 0.014672 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 152 0.277825 0.107643 0.024008 0.000000 0.000000 0.000000 0.000000 0.000000
#> 156 0.376767 0.311824 0.140504 0.037476 0.000000 0.000000 0.000000 0.000000
#> 160 0.205988 0.358291 0.338702 0.176342 0.056012 0.000000 0.000000 0.000000
#> 164 0.047010 0.169985 0.332235 0.356788 0.213350 0.080430 0.000000 0.000000
#> 168 0.004615 0.034404 0.137636 0.300875 0.365191 0.249410 0.111311 0.000000
#> 172 0.000200 0.003057 0.024833 0.109475 0.266492 0.363841 0.282318 0.148903
#> 176 0.000004 0.000123 0.002006 0.017694 0.085631 0.231162 0.353416 0.310035
#> 180 0.000000 0.000002 0.000074 0.001304 0.012455 0.065934 0.196616 0.335187
#> 184 0.000000 0.000000 0.000001 0.000045 0.000841 0.008668 0.050021 0.164165
#> 188 0.000000 0.000000 0.000000 0.000001 0.000027 0.000538 0.005967 0.037421
#> 192 0.000000 0.000000 0.000000 0.000000 0.000000 0.000016 0.000342 0.004067
#> 196 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000010 0.000216
#> 200 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000006
#>          176      180      184      188      192      196 200
#> 80  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 84  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 88  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 92  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 96  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 100 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 104 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 108 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 112 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 116 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 120 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 124 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 128 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 132 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 136 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 140 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 144 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 148 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 152 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 156 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 160 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 164 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 168 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 172 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 176 0.193063 0.000000 0.000000 0.000000 0.000000 0.000000   0
#> 180 0.330907 0.243250 0.000000 0.000000 0.000000 0.000000   0
#> 184 0.310823 0.343852 0.298841 0.000000 0.000000 0.000000   0
#> 188 0.134697 0.282193 0.348745 0.361548 0.000000 0.000000   0
#> 192 0.027629 0.108720 0.251386 0.348386 0.445664 0.000000   0
#> 196 0.002745 0.020149 0.086490 0.221669 0.354374 0.604138   0
#> 200 0.000135 0.001837 0.014538 0.068397 0.199962 0.395862   1
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
