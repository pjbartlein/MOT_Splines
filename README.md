# Introduction

# Introduction

Here we show the R code that we used to reproduce Matlab `csaps()`
smoothing spline fit to the Allen Hills “long shapshot” in Shakleton et
al. (2026), and to fit an alternative optimized smoothing spline in R. A
moothing splines is a cubic spline $``f\left(x_{i}\right)``$ that is fit
to\$data $``y_{i}``$ by minimizing a “penalized” residual sum-of-squares
``` math
RSS\left(f,\lambda \right)=\{\textsf{goodness of fit}\}+\lambda\{\textsf{roughness}\}
```

The first term represents the goodness-of-fit of the spline curve to the
data, while the second “penalty” term describes the rigidity or local
variability of the curve. If λ = 0, the fitted spline simply
interpolates between the $`y_{i}`$’s, while as λ → ∞ the penalty term
increases and the optimal fitted curve becomes a least-squares straight
line (Hastie and Tibshirani, 1990, p. 27).

The objective function used by Matlab `csaps()` function is
``` math
RSS\left(f,\mathrm{\lambda }\right)=\left\{\ p\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\left(1-p\right)\int {\left[f^{\prime \prime }\left(x\right)\right]}^{2}\,dx\right\}
```

Where $`p`$ is the user-selected smoothing parameter, and
$`p=1/\left(1+\mathrm{\lambda }\right)`$. When $`p=0`$, a straight-line
fit is obtained (i.e. as $`λ→∞, p→0`$), while if $`λ=0`$, then $`p=1`$,
and the spline $`(f\left(x_{i}\right))`$ simply interpolates between
points.

# Set up

Set path and filenames:

``` r
# paths for input .csv files -- modify as appropriate
datapath <- "~/Projects/globalT/MOT_Splines/data"
```

``` math
RSS\left(f,\lambda \right)=\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\lambda \int {{\left[f^{\prime \prime }\left(x\right)\right]}^{2}dx,}
```
$`RSS\left(f,\lambda \right)=\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\lambda \int {{\left[f^{\prime \prime }\left(x\right)\right]}^{2}dx,}`$
