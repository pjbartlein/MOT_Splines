library(readxl, writexl)
library(Hmisc)
library(matrixStats)
library(npreg)

# bootstraped splines =============================================================================================

# Allen_Hills_long
data_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/MOTs/Shackleton/"
file_name <- "AH_long.xlsx"
xlxs_file <- paste(data_path, file_name, sep = "")
AH_long <- read_excel(xlxs_file)
head(AH_long); tail(AH_long)

# MOT + 0.5C
AH_long$MOT_anomaly <- AH_long$MOT_anomaly + 0.5

# spline from paper
file_name <- "AH_long_spline_paper.xlsx"
xlxs_file <- paste(data_path, file_name, sep = "")
AH_long_spline_paper <- read_excel(xlxs_file)
head(AH_long_spline_paper); tail(AH_long_spline_paper)


# Shackleton et al. lambda
age_max <- max(AH_long$AgeMa)
age_min <- min(AH_long$AgeMa)
delta_age <- age_max - age_min
nsamp <- length(unique(AH_long$AgeMa))
print(delta_age); print(nsamp)

# Shackleton et al. cut-off period
cop <- 0.4 # [Myr]
dx <- delta_age / nsamp # [Myr]
print(delta_age); print(dx)

# from their .m file
lambda1 <- (cop / (2.0 * pi))^4 / dx
print(lambda1)
p <- 1.0 /(1.0 + lambda1)
print(p)

# target pts for predicted spline
x_targ <- seq(age_min, age_max, by = 0.005)
head(x_targ); tail(x_targ)
nfits <- length(x_targ)
nfits

# set up

weights <- 1/(AH_long$MOT_abs_err^2)
ones = rep(1, length(AH_long$AgeMa))

# bootstrap estimates of fitted values
set.seed(1234)

ncurves <- 3
nreps <- 10000

fits <- matrix(0, nrow = nsamp, ncol = ncurves)
residuals <- matrix(0, nrow = nsamp, ncol = ncurves)
fits_bootstrap <- matrix(0, nrow=nfits, ncol = nreps)
df <- as.numeric(rep(0, ncurves))
lambda <- as.numeric(rep(0, ncurves))
rss <- as.numeric(rep(0, ncurves))
tss <- as.numeric(rep(0, ncurves))
sigma2 <- as.numeric(rep(0, ncurves))
AIC <- as.numeric(rep(0, ncurves))
rsq <- as.numeric(rep(0, ncurves))
adj_rsq <- as.numeric(rep(0, ncurves))
F <- as.numeric(rep(0, ncurves))

model_colors <- c("black", "purple", "red")

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT_anomaly", main = "Allen_Hills_long + 0.5 C")

# add spline from paper
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_anomaly + 0.5, lwd=2, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_lower + 0.5, lwd=1, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_upper + 0.5, lwd=1, col="blue")
leg_label <- paste("Spline from Shackleton et al. 2026 (Fig. 1)", sep="")
legend(3.5, -2.75, legend=c(leg_label), col="blue", lwd=2, bty="n")
#dev_print(pdf, "AH_long_Shackleton_spline.pdf")

# Models

# Shackleton et al. spline

spar <- 0.69
curve_01 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones, spar = spar, keep.stuff = TRUE) # Shackleton et al. specs. for lambda
# curve_01 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones, df = 13) # Best guess
y_pred_curve_01 <- predict(curve_01, x_targ, deriv = 0)$y
y_d1_curve_01 <- predict(curve_01, x_targ, deriv = 1)$y
y_d2_curve_01 <- predict(curve_01, x_targ, deriv = 2)$y
df[1] <- curve_01$df; lambda[1] <- curve_01$lambda; fits[, 1] <- curve_01$y
residuals[, 1] <- AH_long$MOT_anomaly - curve_01$y
sigma2[1] <- sum(residuals[, 1]^2)/(nsamp - 1.0)
rss[1] <- sum(residuals[, 1]^2 )
tss[1] <- sum((AH_long$MOT_anomaly - mean(AH_long$MOT_anomaly))^2)
# http://users.stat.umn.edu/~helwig/notes/smooth-spline-notes.html#aic-and-bic
max_log_lik <- (-1.0 / (2.0 * sigma2[1])) * rss[1] - 1.0 * (nsamp / 2.0) * log(sigma2[1]) - 1.0 * (nsamp / 2.0) * log(2 * pi)
AIC[1] <- - 2.0 * max_log_lik + 2.0 * df[1]
AIC[1] <- round(AIC[1], 2)
print(AIC[1])

# R-sq values
rsq[1] <- 1 - (rss[1]/tss[1])
adj_rsq[1] <- 1 - (rss[1] / (nsamp - df[1])) / (tss[1] / (nsamp - 1))
rsq[1] <- round(rsq[1], digits = 3)
adj_rsq[1] <- round(adj_rsq[1], digits = 3)
rsq[1]; adj_rsq[1]

hist(residuals[, 1], breaks = seq(-4.25, 4.25, by = 0.25))
qqnorm(residuals[, 1]); qqline(residuals[, 1])
shapiro.test(residuals[, 1])
acf(residuals[, 1])
Box.test(residuals[, 1], lag = 15, type = "Ljung-Box", fitdf = df[1])

oldpar <- par(pty="s")
ylim_resid <- c(-2, 2)
plot(residuals[, 1] ~ fits[, 1], ylim = ylim_resid, xlim = ylim_resid, pch = 16, cex = 0.8, 
     ylab = "Residuals from Smoothing Spline", xlab =  "Smoothing Spline Fitted Values")
lo_resid <- loess(residuals[, 1] ~ fits[, 1], span=0.80, degree = 1, iterations = 3)
lines(lo_resid$fitted ~ fits[, 1], col = "red", lwd = 2)
par(oldpar)

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

# bootsrap CL's
fits_bootstrap <- matrix(0, nrow=nfits, ncol = nreps)
k <- 0; j <- 1
for (i in 1:nreps) {
  # print(i)
  x_boot <- AH_long$AgeMa + rnorm(nsamp, mean = 0, sd = AH_long$Rel_err_Myr)
  x_boot <- x_boot[order(x_boot)]
  y_boot <- AH_long$MOT_anomaly + rnorm(nsamp, mean = 0, sd = AH_long$MOT_abs_err)
  y_boot <- y_boot[order(x_boot)]
  # points(x_boot, y_boot, pch = 16, cex = 0.5, col = "red")
  
  fit_boot <- smooth.spline(x_boot, y_boot, w = ones, spar = spar ) # df = df[1])
  if (length(x_boot) == length(fit_boot$y)) {
    y_pred <- predict(fit_boot, x_targ, deriv = 0)$y
    if (i < 500) {
      lines(x_targ, y_pred, col=adjustcolor(model_colors[j], alpha = 0.05))
      # points(AH_long$AgeMa, y_pred, pch = 16, cex = 0.5, col = "red")
    }
    k <- k + 1
    fits_bootstrap[1:nfits, k] <- y_pred[1:nfits]
  }
}

boot_cl_01 <- rowQuantiles(fits_bootstrap, probs = c(0.025, 0.16, 0.50, 0.84, 0.975))
head(boot_cl_01)
lines(x_targ, boot_cl_01[, 2], col = model_colors[1])
lines(x_targ, boot_cl_01[, 4], col = model_colors[1])
lines(x_targ, boot_cl_01[, 1], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 5], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 3], lty = "dashed", col = model_colors[1])

# Shackleton et al. specs. for lambda
lines(x_targ, y_pred_curve_01, lwd=2, col=model_colors[1])
leg_label <- paste("Smoothing Spline, spar = ", as.character(round(spar, 3)), ", lambda = ", as.character(round(lambda[1],7)), 
                   ", df = ", as.character(round(df[1], 2)), ", AIC = ", as.character(round(AIC[1], 2)),sep="")
legend(3.5, -2.75, legend=c(leg_label), col=c(model_colors[1]), lwd=2, bty="n")

#dev_print(pdf, "AH_long_Shackleton_spline.pdf")

# spline from paper

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

lines(x_targ, boot_cl_01[, 2], col = model_colors[1])
lines(x_targ, boot_cl_01[, 4], col = model_colors[1])
lines(x_targ, boot_cl_01[, 1], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 5], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 3], lty = "dashed", col = model_colors[1])
lines(x_targ, y_pred_curve_01, lwd=2, col=model_colors[1])

lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_anomaly + 0.5, lwd=2, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_lower + 0.5, lwd=1, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_upper + 0.5, lwd=1, col="blue")
leg_label <- paste("Spline from Shackleton et al. 2026 (Fig. 1)", sep="")
legend(3.5, -2.50, legend=c(leg_label), col="blue", lwd=2, bty="n")

points(AH_long$MOT_anomaly ~ AH_long$AgeMa, pch = 16, cex = 0.8, col = "blue")

#dev_print(pdf, "AH_long_spline_comparison.pdf")

# fitted values

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "gray",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

lines(fits[, 1] ~ AH_long$AgeMa, lwd = 0.5, col = "black")
points(fits[, 1] ~ AH_long$AgeMa, pch = 16, cex = 0.4, col = "black")
points(y_pred_curve_01 ~ x_targ, pch = "+", cex = 0.15, col = "purple")

#dev_print(pdf, "AH_long_fitted_comparison.pdf")

# derivatives
plot(y_d1_curve_01 ~ x_targ, type = "l", lwd = 1, col = "red")
plot(y_d2_curve_01 ~ x_targ, type = "l", lwd = 1, col = "red")

# optimal

curve_02 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones) # estimate lambda, spar, df
df[2] <- curve_02$df; lambda[2] <- curve_02$lambda; fits[, 2] <- curve_02$y
y_pred_curve_02 <- predict(curve_02, x_targ, deriv = 0)$y
y_d1_curve_02 <- predict(curve_02, x_targ, deriv = 1)$y
y_d2_curve_02 <- predict(curve_02, x_targ, deriv = 2)$y
residuals[, 2] <- AH_long$MOT_anomaly - curve_02$y
sigma2[2] <- sum(residuals[, 2]^2)/(nsamp - 1.0)
rss[2] <- sum(residuals[, 2]^2 )
tss[2] <- sum((AH_long$MOT_anomaly - mean(AH_long$MOT_anomaly))^2)
# http://users.stat.umn.edu/~helwig/notes/smooth-spline-notes.html#aic-and-bic
max_log_lik <- (-1.0 / (2.0 * sigma2[2])) * rss[2] - 1.0 * (nsamp / 2.0) * log(sigma2[2]) - 1.0 * (nsamp / 2.0) * log(2 * pi)
AIC[2] <- - 2.0 * max_log_lik + 2.0 * df[2]
AIC[2] <- round(AIC[2], 2)
print(AIC[2])

# R-sq values
rsq[2] <- 1 - (rss[2]/tss[2])
adj_rsq[2] <- 1 - (rss[1] / (nsamp - df[2])) / (tss[2] / (nsamp - 1))
rsq[2] <- round(rsq[2], digits = 3)
adj_rsq[2] <- round(adj_rsq[2], digits = 3)
rsq[2]; adj_rsq[2]

hist(residuals[, 2], breaks = seq(-4.25, 4.25, by = 0.25))
qqnorm(residuals[, 2]); qqline(residuals[, 2])
shapiro.test(residuals[, 2])
acf(residuals[, 2])
Box.test(residuals[, 2], lag = 15, type = "Ljung-Box", fitdf = df[2])

oldpar <- par(pty="s")
ylim_resid <- c(-2, 2)
plot(residuals[, 2] ~ fits[, 2], ylim = ylim_resid, xlim = ylim_resid, pch = 16, cex = 0.8, 
    ylab = "Residuals from Smoothing Spline", xlab =  "Smoothing Spline Fitted Values")
lo_resid <- loess(residuals[, 2] ~ fits[, 2], span=0.80, degree = 1, iterations = 3)
lines(lo_resid$fitted ~ fits[, 2], col = "red", lwd = 2)
par(oldpar)

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT_anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

# bootsrap CL's
fits_bootstrap <- matrix(0, nrow=nfits, ncol = nreps)
k <- 0; j <- 2
for (i in 1:nreps) {
  # print(i)
  x_boot <- AH_long$AgeMa + rnorm(nsamp, mean = 0, sd = AH_long$Rel_err_Myr)
  x_boot <- x_boot[order(x_boot)]
  y_boot <- AH_long$MOT_anomaly + rnorm(nsamp, mean = 0, sd = AH_long$MOT_abs_err)
  y_boot <- y_boot[order(x_boot)]
  # points(x_boot, y_boot, pch = 16, cex = 0.5, col = "red")
  
  fit_boot <- smooth.spline(x_boot, y_boot, w = ones, df = df[2])
  if (length(x_boot) == length(fit_boot$y)) {
    y_pred <- predict(fit_boot, x_targ, deriv = 0)$y
    if (i < 500) {
      lines(x_targ, y_pred, col=adjustcolor(model_colors[2], alpha = 0.02))
      # points(AH_long$AgeMa, y_pred, pch = 16, cex = 0.5, col = "red")
    }
    k <- k + 1
    fits_bootstrap[1:nfits, k] <- y_pred[1:nfits]
  }
}

boot_cl_02 <- rowQuantiles(fits_bootstrap, probs = c(0.025, 0.16, 0.50, 0.84, 0.975))
head(boot_cl_01)
lines(x_targ, boot_cl_02[, 2], col = model_colors[2])
lines(x_targ, boot_cl_02[, 4], col = model_colors[2])
lines(x_targ, boot_cl_02[, 1], lty = "dotted", lwd = 2, col = model_colors[2])
lines(x_targ, boot_cl_02[, 5], lty = "dotted", lwd = 2, col = model_colors[2])
lines(x_targ, boot_cl_02[, 3], lty = "dashed", col = model_colors[2])

# optimal lambda
lines(x_targ, y_pred_curve_02, lwd=2, col=model_colors[2])
leg_label <- paste("Smoothing Spline, estimated lambda = ", as.character(round(lambda[2],7)), 
                   ", df = ", as.character(round(df[2], 2)), ", AIC = ", as.character(round(AIC[2], 2)),sep="")
legend(3.5, -2.75, legend=c(leg_label), col=c(model_colors[2]), lwd=2, bty="n")

points(AH_long$MOT_anomaly ~ AH_long$AgeMa, pch = 16, cex = 0.8, col = "blue")
#dev_print(pdf, "AH_long_estimated_spline.pdf")

# derivatives
plot(y_d1_curve_02 ~ x_targ, type = "l", lwd = 1, col = "red")
plot(y_d2_curve_02 ~ x_targ, type = "l", lwd = 1, col = "red")

# OLS

curve_03 <- lm(AH_long$MOT_anomaly ~ AH_long$AgeMa) # OLS
summary(curve_03)
# y_pred_curve_03 <- predict(curve_03, data.frame(x_targ))
y_pred_curve_03 <- curve_03$coefficients[1] + curve_03$coefficients[2] * x_targ
fits[, 3] <- curve_03$fitted.values
residuals[, 3] <- AH_long$MOT_anomaly - curve_03$fitted.values
rss[3] <- sum(residuals[, 3]^2 )
tss[3] <- sum((AH_long$MOT_anomaly - mean(AH_long$MOT_anomaly))^2)
AIC[3] <- round(AIC(curve_03), 2)
print(AIC[3])

# R-sq values
rsq[3] <- 1 - (rss[3]/tss[3])
adj_rsq[3] <- 1 - (rss[3] / (nsamp - df[3])) / (tss[3] / (nsamp - 1))
rsq[3] <- round(rsq[3], digits = 3)
adj_rsq[3] <- round(adj_rsq[3], digits = 3)
rsq[3]; adj_rsq[3]

hist(residuals[, 3], breaks = seq(-4.25, 4.25, by = 0.25))
qqnorm(residuals[, 3]); qqline(residuals[, 3])
shapiro.test(residuals[, 3])
acf(residuals[, 3])
Box.test(residuals[, 3], lag = 15, type = "Ljung-Box", fitdf = df[3])

oldpar <- par(pty="s")
ylim_resid <- c(-2, 2)
plot(residuals[, 3] ~ fits[, 3], ylim = ylim_resid, xlim = ylim_resid, pch = 16, cex = 0.8, 
     ylab = "Residuals from Smoothing Spline", xlab =  "Smoothing Spline Fitted Values")
lo_resid <- loess(residuals[, 3] ~ fits[, 3], span=0.80, degree = 1, iterations = 3)
lines(lo_resid$fitted ~ fits[, 3], col = "red", lwd = 2)
par(oldpar)

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

# bootsrap CL's
fits_bootstrap <- matrix(0, nrow=nfits, ncol = nreps)
k <- 0; j <- 3
for (i in 1:nreps) {
  # print(i)
  x_boot <- AH_long$AgeMa + rnorm(nsamp, mean = 0, sd = AH_long$Rel_err_Myr)
  x_boot <- x_boot[order(x_boot)]
  y_boot <- AH_long$MOT_anomaly + rnorm(nsamp, mean = 0, sd = AH_long$MOT_abs_err)
  y_boot <- y_boot[order(x_boot)]
  # points(x_boot, y_boot, pch = 16, cex = 0.5, col = "red")
  
  fit_boot <- lm(y_boot ~ x_boot) # OLS
  if (length(x_boot) == length(fit_boot$fitted.values)) {
    y_pred <- fit_boot$coefficients[1] + fit_boot$coefficients[2] * x_targ
    if (i < 500) {
      lines(x_targ, y_pred, col=adjustcolor(model_colors[3], alpha = 0.02))
      # points(x_targ, y_pred, pch = 16, cex = 0.5, col = "red")
    }
    k <- k + 1
    fits_bootstrap[1:nfits, k] <- y_pred[1:nfits]
  }
}
boot_cl_03 <- rowQuantiles(fits_bootstrap, probs = c(0.025, 0.16, 0.50, 0.84, 0.975))
head(boot_cl_03)
lines(x_targ, boot_cl_03[, 2], col = model_colors[3])
lines(x_targ, boot_cl_03[, 4], col = model_colors[3])
lines(x_targ, boot_cl_03[, 1], lty = "dotted", lwd = 2, col = model_colors[3])
lines(x_targ, boot_cl_03[, 5], lty = "dotted", lwd = 2, col = model_colors[3])
lines(x_targ, boot_cl_03[, 3], lty = "dashed", col = model_colors[3])

lines(x_targ, y_pred_curve_03, lwd=2, col=model_colors[3])
leg_label <- paste("OLS Linear Model", ", AIC = ", as.character(round(AIC[3], 2)), sep="")
legend(3.5, -2.75, legend=c(leg_label), col=model_colors[3], lwd=2, bty="n")

points(AH_long$MOT_anomaly ~ AH_long$AgeMa, pch = 16, cex = 0.8, col = "blue")
#dev_print(pdf, "AH_long_OLS.pdf")

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

lines(x_targ, boot_cl_01[, 2], col = model_colors[1])
lines(x_targ, boot_cl_01[, 4], col = model_colors[1])
lines(x_targ, boot_cl_01[, 1], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 5], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, y_pred_curve_01, lwd=2, col=model_colors[1])
leg_label <- paste("Smoothing Spline, lambda = ", as.character(round(lambda[1],7)), 
                   ", df = ", as.character(round(df[1], 2)), ", AIC = ", as.character(round(AIC[1], 2)),sep="")
legend(3.5, -2.25, legend=c(leg_label), col=c(model_colors[1]), lwd=2, bty="n")

lines(x_targ, boot_cl_02[, 2], col = model_colors[2])
lines(x_targ, boot_cl_02[, 4], col = model_colors[2])
lines(x_targ, boot_cl_02[, 1], lty = "dotted", lwd = 2, col = model_colors[2])
lines(x_targ, boot_cl_02[, 5], lty = "dotted", lwd = 2, col = model_colors[2])
lines(x_targ, y_pred_curve_02, lwd=2, col=model_colors[2])
leg_label <- paste("Smoothing Spline, estimated lambda = ", as.character(round(lambda[2],7)), 
                   ", df = ", as.character(round(df[2], 2)), ", AIC = ", as.character(round(AIC[2], 2)),sep="")
legend(3.5, -2.50, legend=c(leg_label), col=c(model_colors[2]), lwd=2, bty="n")

lines(x_targ, boot_cl_03[, 2], col = model_colors[3])
lines(x_targ, boot_cl_03[, 4], col = model_colors[3])
lines(x_targ, boot_cl_03[, 1], lty = "dotted", lwd = 2, col = model_colors[3])
lines(x_targ, boot_cl_03[, 5], lty = "dotted", lwd = 2, col = model_colors[3])
lines(x_targ, y_pred_curve_03, lwd=2, col=model_colors[3])
leg_label <- paste("OLS Linear Model", ", AIC = ", as.character(round(AIC[3], 2)), sep="")
legend(3.5, -2.75, legend=c(leg_label), col=model_colors[3], lwd=2, bty="n")
#dev_print(pdf, "AH_long_all.pdf")
#dev_print(pdf, "AH_long_est_vs_OLD.pdf")

# bootstrap uncertainties
plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

lines(AH_long$AgeMa, fits[, 1], lwd = 2, col = model_colors[1])
lines(AH_long$AgeMa, fits[, 2], lwd = 2, col = model_colors[2])
lines(AH_long$AgeMa, fits[, 3], lwd = 2, col = model_colors[3])

out_data_uncert <- data.frame(x_targ, 
                       y_pred_curve_01, boot_cl_01[, 1], boot_cl_01[, 2], boot_cl_01[, 3], boot_cl_01[, 4], boot_cl_01[, 5],
                       y_pred_curve_02, boot_cl_02[, 1], boot_cl_02[, 2], boot_cl_02[, 3], boot_cl_02[, 4], boot_cl_02[, 5],
                       y_pred_curve_03, boot_cl_03[, 1], boot_cl_03[, 2], boot_cl_03[, 3], boot_cl_03[, 4], boot_cl_03[, 5] )
head(out_data_uncert)
names(out_data_uncert) <- c("AgeMA", 
                     "Shack_fit", "Shack_q0.025", "Shack_q0.16", "Shack_median", "Shack_q0.84", "Shack_q0.975",
                     "opt_fit", "opt_q0.025", "opt_q0.16", "opt_median", "opt_q0.84", "opt_q0.975",
                     "linear_fit", "linear_q0.025", "linear_q0.16", "linear_median", "linear_q0.84", "linear_q0.975" )
head(out_data_uncert)
write.csv(out_data_uncert, "AH+0.5_splines_uncert.csv", row.names = FALSE)

# fitted at target points
plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

lines(x_targ, y_pred_curve_01, lwd = 2, col = model_colors[1])
lines(x_targ, y_pred_curve_02, lwd = 2, col = model_colors[2])
lines(x_targ, y_pred_curve_03, lwd = 2, col = model_colors[3])

out_data_fitted <- data.frame(x_targ, y_pred_curve_01, y_pred_curve_02, y_pred_curve_03 )
head(out_data_fitted)
names(out_data_fitted) <- c("AgeMA", "Shackleton_fitted", "optimized_spline", "OLS")
head(out_data_fitted)
write.csv(out_data_fitted, "AH+0.5_splines_fitted_at_target_points.csv", row.names = FALSE)

# fitted at data points
plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

lines(AH_long$AgeMa, fits[, 1], lwd = 2, col = model_colors[1])
lines(AH_long$AgeMa, fits[, 2], lwd = 2, col = model_colors[2])
lines(AH_long$AgeMa, fits[, 3], lwd = 2, col = model_colors[3])

out_data_fitted <- data.frame(AH_long$AgeMa, fits[, 1], fits[, 2], fits[, 3] )
head(out_data_fitted)
names(out_data_fitted) <- c("AgeMA", "Shackleton_fitted", "optimized_spline", "OLS")
head(out_data_fitted)
write.csv(out_data_fitted, "AH+0.5_splines_fitted_at_data_points.csv", row.names = FALSE)



# compare models
print(AIC[1:3])
print(rsq[1:3])

# F statistics
F[1] <- ((tss[1] - rss[1])/(df[1] - 1)) / (rss[1]/(nsamp - df[1]))
F[2] <- ((tss[2] - rss[2])/(df[2] - 1)) / (rss[2]/(nsamp - df[2]))
F[3] <- summary(curve_03)$fstatistic[1]
F <- as.numeric(F)
print(F)

# p-values
p1 <- pf(F[1], df[1] - 1, nsamp - df[1], lower.tail = FALSE)
p2 <- pf(F[2], df[2] - 1, nsamp - df[2], lower.tail = FALSE)
p3 <- pf(F[3], 2, nsamp - 1, lower.tail = FALSE)
print(c(p1, p2, p3))

# compare 1 and 2
F_12 <- ((rss[1] - rss[2]) / (df[2] - df[1])) / (rss[2] / (nsamp - df[2]))
# F_12 <- ((rss[2] - rss[1]) / (df[1] - df[2])) / (rss[1] / (nsamp - df[1]))
p12 <- pf(F_12, (df[1] - df[2]), (nsamp - df[1]), lower.tail = FALSE)

library(mgcv)
        
curve_04 <- gam(AH_long$MOT_anomaly ~ s(AH_long$AgeMa), method = "REML")
summary(curve_04)


opar <- par(mar=c(5,4,4,5)+0.1) # space for second axis

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-4.0, 2.0), xlab = "Age Ma", ylab = "MOT anomaly", 
     font.main = 1, main = "Allen Hills (long) + 0.5 C")
minor.tick(nx = 10, tick.ratio = 0.5)

lines(x_targ, boot_cl_01[, 2], col = model_colors[1])
lines(x_targ, boot_cl_01[, 4], col = model_colors[1])
# lines(x_targ, boot_cl_01[, 1], lty = "dotted", lwd = 2, col = model_colors[1])
# lines(x_targ, boot_cl_01[, 5], lty = "dotted", lwd = 2, col = model_colors[1])
# lines(x_targ, boot_cl_01[, 3], lty = "dashed", col = model_colors[1])
lines(x_targ, y_pred_curve_01, lwd=2, col=model_colors[1])

par(new=TRUE)                   # second plot is going to get added to first
# plot(spline_log10dx$x, log10(spline_log10dx$y), lwd=2, type = "l", col="palegreen3", axes = FALSE,
#      xlim = c(3.5, 0), ylim = c(1., 3.), xlab = "Age Ma", ylab = "MOT anomaly")
# plot(x_targ, log10(y_spline_log10dx_2), lwd=2, type = "l", col="palegreen3", axes = FALSE,
#      xlim = c(3.5, 0), ylim = c(1., 3.), xlab = "Age Ma", ylab = "MOT anomaly")
plot(log10(AH_long2$dx) ~ AH_long2$AgeMa,  cex = 0.8, pch = 16, col = "green",
     xlim = c(3.5, 0), ylim = c(1., 3.), xlab = " ", ylab = " ", axes = FALSE)
lines(x_targ, log10(y_spline_log10dx_2), lwd = 2, col = "palegreen3")
axis(side=4)                    # add axis
mtext(side=4,line=3.8,"log10(dx)")   # add label
par(opar)

#dev_print(pdf, "AH_long_dx_vs_fit.pdf")


# y_curve_01 <- predict(curve_01, AH_long2$AgeMa, deriv = 0)$y
# y_spline_log10dx <- predict(spline_log10dx, AH_long2$AgeMa, deriv = 0)$y
# plot(y_spline_log10dx, y_curve_01)
# cor(y_spline_log10dx, y_curve_01)
# 
# y_spline_log10dx_2 <- predict(spline_log10dx, x_targ, deriv = 0)$y
# plot(y_spline_log10dx_2, y_pred_curve_01)
# cor(y_spline_log10dx_2, y_pred_curve_01)
# 
# lm_ypred <- lm(y_pred_curve_01 ~ x_targ)
# lm_dx <- lm(y_spline_log10dx_2 ~ x_targ)
# 
# plot(lm_dx$residuals, lm_ypred$residuals)
