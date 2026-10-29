
Here we show the R code that we used to reproduce Matlab `csaps()` smoothing spline fit to the Allen Hills “long shapshot” in Shackleton et
al. (2026), and to fit an alternative optimized smoothing spline in R. A smoothing splines is a cubic spline $f\left(x_{i}\right)$ that is fit
to some data $y_{i}$ by minimizing a “penalized” residual sum-of-squares

$$
RSS\left(f,\lambda \right)=\{\textsf{goodness of fit}\}+\lambda\{\textsf{roughness}\}
$$

The first term represents the goodness-of-fit of the spline curve to the data, while the second “penalty” term describes the rigidity or local
variability of the curve. If λ = 0, the fitted spline simply interpolates between the $y_{i}$’s, while as λ → ∞ the penalty term
increases and the optimal fitted curve becomes a least-squares straight line (Hastie and Tibshirani, 1990, p. 27).

The objective function used by the Matlab `csaps()` function is
$$
RSS\left(f,\mathrm{\lambda }\right)=\left\{\ p\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\left(1-p\right)\int {\left[f^{\prime \prime }\left(x\right)\right]}^{2}\,dx\right\}
$$

Where $f^{\prime \prime }\left(x\right)$ is the second derivative of the spline curve, and the $w_{i}$’s are possible weights for each data value, $p$ is the user-selected smoothing parameter, and $p=1/\left(1+\mathrm{\lambda }\right)$. When $p=0$, a straight-line
fit is obtained (i.e. as $λ→∞$, $p→0$), while if $λ=0$, then $p=1$ and the spline $(f\left(x_{i}\right))$ simply interpolates between
points.

The objective function for the R `smooth.spline()` function is
$$
RSS\left(f,\lambda \right)=\sum_{i=1}^{n} {w_{i}{\left(y_{i}-f\left(x_{i}\right)\right)}^{2}}+\lambda \int {{\left[f^{\prime \prime }\left(x\right)\right]}^{2}dx,}
$$

Again, If $λ=0$, the fitted spline simply interpolates between the $y_{i}$’s, while as $λ→∞$ the penalty term increases and the optimal fitted curve becomes a least-squares straight line. The smoothing parameter is thus a measure of the “wiggliness” of the curve and can be set by a user at some fixed, scientifically based value, or it can be estimated by Generalized Cross Validation (GCV, Wahba, 1990, sec. 4.3) or some other approach. (GCV is a computationally efficient “leave-one-out” approach that is used to find an optimum trade-off between the goodness-of-fit and wigglyness of smoothing spline: an extremely smooth spline won’t fit the data very well, while an extremely wiggly spline can fit the data perfectly at the expense of being extremely non-smooth.) We use this approach below to obtain and "optimal" smoothing spline.

In the case of the R function, there are three ways of specifying or describing the amount of smoothing: via $λ$; via a *spar* value, a scaleless description of the amount of smoothing which generally varies between 0 an 1; and via the equivalent degrees of freedom, $df$ (Hastie and Tibshirani, 1990, sec. 3.5), a general measure of the amount of fitting done by the spline. Low values of $df$ describe a fitted curve with fewer inflections while higher values describe a more wiggly curve.

Because the internal algorithms in the R and Matlab functions differ, there is no direct conversion between the methods (i.e. selecting the same value of $λ$ for both will not yield equivalent results). The R function can nevertheless be used to approximate the results of applying the Matlab function by iteratively varying the *spar* parameter, and minimizing the sum-of-squared differences between the R- and Matlab-fitted smoothing splines. We use this approach below to approximate the results of from the Matlab `csaps()` using the R `smooth.spline()` functio
