# Introduction

Here we show the `R` code that we used to reproduce Matlab `csaps()`
smoothing spline fit to the Allen Hills “long shapshot” in Shakleton et
al. (2026).

# Set up

Set path and filenames:

``` r
# paths for input .csv files -- modify as appropriate
datapath <- "~/Projects/globalT/MOT_Splines/data"
```

$$RSS\left(f,\lambda \right)=\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\lambda \int {{\left[f^{\prime \prime }\left(x\right)\right]}^{2}dx,}$$

``` math
RSS\left(f,\lambda \right)=\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\lambda \int {{\left[f^{\prime \prime }\left(x\right)\right]}^{2}dx,}
```
$`RSS\left(f,\lambda \right)=\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\lambda \int {{\left[f^{\prime \prime }\left(x\right)\right]}^{2}dx,}`$
