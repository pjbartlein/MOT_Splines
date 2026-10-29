# Introduction

Here we show the R code that we used to reproduce Matlab `csaps()`
smoothing spline fit to the Allen Hills “long shapshot” in Shackleton et
al. (2026), and to fit an alternative optimized smoothing spline in R. A
smoothing splines is a cubic spline $`f\left(x_{i}\right)`$ that is fit
to some data $`y_{i}`$ by minimizing a “penalized” residual
sum-of-squares

``` math
RSS\left(f,\lambda \right)=\{\textsf{goodness of fit}\}+\lambda\{\textsf{roughness}\}
```

The first term represents the goodness-of-fit of the spline curve to the
data, while the second “penalty” term describes the rigidity or local
variability of the curve. If λ = 0, the fitted spline simply
interpolates between the $`y_{i}`$’s, while as λ → ∞ the penalty term
increases and the optimal fitted curve becomes a least-squares straight
line (Hastie and Tibshirani, 1990, p. 27).

The objective function used by the Matlab `csaps()` function is
``` math
RSS\left(f,\mathrm{\lambda }\right)=\left\{\ p\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\left(1-p\right)\int {\left[f^{\prime \prime }\left(x\right)\right]}^{2}\,dx\right\}
```

Where $`f^{\prime \prime }\left(x\right)`$ is the second derivative of
the spline curve, and the $`w_{i}`$’s are possible weights for each data
value, $`p`$ is the user-selected smoothing parameter, and
$`p=1/\left(1+\mathrm{\lambda }\right)`$. When $`p=0`$, a straight-line
fit is obtained (i.e. as $`λ→∞`$, $`p→0`$), while if $`λ=0`$, then
$`p=1`$ and the spline $`(f\left(x_{i}\right))`$ simply interpolates
between points.

The objective function for the R `smooth.spline()` function is
``` math
RSS\left(f,\lambda \right)=\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\lambda \int {{\left[f^{\prime \prime }\left(x\right)\right]}^{2}dx,}
```

Again, If $`λ=0`$, the fitted spline simply interpolates between the
$`y_{i}`$’s, while as $`λ→∞`$ the penalty term increases and the optimal
fitted curve becomes a least-squares straight line. The smoothing
parameter is thus a measure of the “wiggliness” of the curve and can be
set by a user at some fixed, scientifically based value, or it can be
estimated by Generalized Cross Validation (GCV, Wahba, 1990, sec. 4.3)
or some other approach. (GCV is a computationally efficient
“leave-one-out” approach that is used to find an optimum trade-off
between the goodness-of-fit and wigglyness of smoothing spline: an
extremely smooth spline won’t fit the data very well, while an extremely
wiggly spline can fit the data perfectly at the expense of being
extremely non-smooth.) We use this approach below to obtain and
“optimal” smoothing spline.

In the case of the R function, there are three ways of specifying or
describing the amount of smoothing: via $`λ`$; via a *spar* value, a
scaleless description of the amount of smoothing which generally varies
between 0 an 1; and via the equivalent degrees of freedom, $`df`$
(Hastie and Tibshirani, 1990, sec. 3.5), a general measure of the amount
of fitting done by the spline. Low values of $`df`$ describe a fitted
curve with fewer inflections while higher values describe a more wiggly
curve.

Because the internal algorithms in the R and Matlab functions differ,
there is no direct conversion between the methods (i.e. selecting the
same value of $`λ`$ for both will not yield equivalent results). The R
function can nevertheless be used to approximate the results of applying
the Matlab function by iteratively varying the *spar* parameter, and
minimizing the sum-of-squared differences between the R- and
Matlab-fitted smoothing splines. We use this approach below to
approximate the results of from the Matlab `csaps()` using the R
`smooth.spline()` functio

# Read data

Load the neccessary packages. `{readxl}` and `{writexl}` read and write
Excel files, `{Hmisc}` adds supplemental tick marks and plots, and
`{matrixStats}` is used to get the bootstrap confidence itervals.

``` r
library(readxl, writexl)
library(Hmisc)
```

    ## 
    ## Attaching package: 'Hmisc'

    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units

``` r
library(matrixStats)
```

Read the Allen Hills “long” snapshot data, and add 0.5 °C to account for
the bias in the MOT anomalies.

``` r
# Allen_Hills_long data
# columns Q, X, Y and Z from the "Data Summary" tab of AllanHills_XeKr_temp-recon_2026-pubData-2.xlsx
data_path <- "~/Projects/globalT/MOT_Splines/data/"
file_name <- "AH_long.xlsx" 
xlxs_file <- paste(data_path, file_name, sep = "")
AH_long <- read_excel(xlxs_file)
head(AH_long); tail(AH_long)
```

    ## # A tibble: 6 × 9
    ##   Age_ka AgeMa Rel_err_Myr Abs_err_Myr MOT_anomaly MOT_rel_err MOT_abs_err Snapshot
    ##    <dbl> <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <chr>   
    ## 1   430  0.43       0.0909      0.102       -1.99        0.559       0.750 'long'  
    ## 2   499. 0.499      0.104       0.118       -0.870       0.431       0.660 'long'  
    ## 3   500  0.5        0.0909      0.106       -2.09        0.574       0.761 'long'  
    ## 4   552. 0.552      0.0618      0.0866      -1.89        0.526       0.726 'long'  
    ## 5   558. 0.558      0.0618      0.0871      -2.04        0.547       0.741 'long'  
    ## 6   579. 0.579      0.0618      0.0888      -1.85        0.519       0.721 'long'  
    ## # ℹ 1 more variable: seqnum <dbl>

    ## # A tibble: 6 × 9
    ##   Age_ka AgeMa Rel_err_Myr Abs_err_Myr MOT_anomaly MOT_rel_err MOT_abs_err Snapshot
    ##    <dbl> <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl> <chr>   
    ## 1  2690   2.69      0.0909       0.310      -0.31        0.313       0.590 'long'  
    ## 2  2848.  2.85      0.0532       0.318       0.495       0.268       0.567 'long'  
    ## 3  2868.  2.87      0.0593       0.321       0.229       0.247       0.558 'long'  
    ## 4  2891.  2.89      0.0593       0.323       1.03        0.319       0.593 'long'  
    ## 5  2895.  2.90      0.0593       0.324       1.38        0.367       0.620 'long'  
    ## 6  3065.  3.06      0.0618       0.343       0.229       0.275       0.570 'long'  
    ## # ℹ 1 more variable: seqnum <dbl>

``` r
# MOT_anomaly + 0.5C
AH_long$MOT_anomaly <- AH_long$MOT_anomaly + 0.5
```

Read Shackleton et al.’s fitted spline, and add 0.5 °C.

``` r
# spline from paper
# columns A-D from the "Splined Records" tab of AllanHills_XeKr_temp-recon_2026-pubData-2.xlsx
file_name <- "AH_long_spline_paper.xlsx" 
xlxs_file <- paste(data_path, file_name, sep = "")
AH_long_spline_paper <- read_excel(xlxs_file)
head(AH_long_spline_paper); tail(AH_long_spline_paper)
```

    ## # A tibble: 6 × 5
    ##   Age_ka AgeMa MOT_anomaly MOT_lower MOT_upper
    ##    <dbl> <dbl>       <dbl>     <dbl>     <dbl>
    ## 1    500 0.5         -1.79     -2.20     -1.39
    ## 2    505 0.505       -1.78     -2.19     -1.39
    ## 3    510 0.51        -1.78     -2.18     -1.40
    ## 4    515 0.515       -1.78     -2.17     -1.40
    ## 5    520 0.52        -1.78     -2.16     -1.41
    ## 6    525 0.525       -1.77     -2.15     -1.41

    ## # A tibble: 6 × 5
    ##   Age_ka AgeMa MOT_anomaly MOT_lower MOT_upper
    ##    <dbl> <dbl>       <dbl>     <dbl>     <dbl>
    ## 1   2975  2.98       0.730     0.316      1.28
    ## 2   2980  2.98       0.721     0.306      1.29
    ## 3   2985  2.98       0.717     0.297      1.30
    ## 4   2990  2.99       0.710     0.286      1.31
    ## 5   2995  3.00       0.701     0.276      1.32
    ## 6   3000  3          0.691     0.263      1.34

``` r
# MOT + 0.5C
AH_long_spline_paper$MOT_anomaly <- AH_long_spline_paper$MOT_anomaly + 0.5
AH_long_spline_paper$MOT_lower <- AH_long_spline_paper$MOT_lower + 0.5
AH_long_spline_paper$MOT_upper <- AH_long_spline_paper$MOT_upper + 0.5
```

Plot the data, and the spline from the paper.

``` r
# Plot the data and spline from the paper
# data
plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age (Ma)", ylab = "MOT anomaly (°C)", 
     font.main = 1, main = "Allen Hills (long) + 0.5 °C")
minor.tick(nx = 10, tick.ratio = 0.5)

# add spline from paper
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_anomaly, lwd=2, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_lower, lwd=1, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_upper, lwd=1, col="blue")
leg_label <- paste("Spline from Shackleton et al. 2026 (Fig. 1)", sep="")
legend(3.5, -2.75, legend=c(leg_label), col="blue", lwd=2, bty="n")
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#dev.print(pdf, "AH_long_Shackleton_spline.pdf")
```
