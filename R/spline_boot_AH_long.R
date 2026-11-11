library(readxl, writexl)
library(Hmisc)
library(matrixStats)

# bootstraped splines =============================================================================================

# Allen_Hills_long data
# columns Q, X, Y and Z from the "Data Summary" tab of AllanHills_XeKr_temp-recon_2026-pubData-2.xlsx
data_path <- "~/Projects/globalT/MOT_Splines/data/"
file_name <- "AH_long.xlsx" 
xlxs_file <- paste(data_path, file_name, sep = "")
AH_long <- read_excel(xlxs_file)
head(AH_long); tail(AH_long)

# MOT_anomaly + 0.5C
AH_long$MOT_anomaly <- AH_long$MOT_anomaly + 0.5

# spline from paper
# columns A-D from the "Splined Records" tab of AllanHills_XeKr_temp-recon_2026-pubData-2.xlsx
file_name <- "AH_long_spline_paper.xlsx" 
xlxs_file <- paste(data_path, file_name, sep = "")
AH_long_spline_paper <- read_excel(xlxs_file)
head(AH_long_spline_paper); tail(AH_long_spline_paper)

# MOT + 0.5C
AH_long_spline_paper$MOT_anomaly <- AH_long_spline_paper$MOT_anomaly + 0.5
AH_long_spline_paper$MOT_lower <- AH_long_spline_paper$MOT_lower + 0.5
AH_long_spline_paper$MOT_upper <- AH_long_spline_paper$MOT_upper + 0.5

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
#dev.print(pdf, "AH_long_Shackleton_spline.pdf")

# Set-up tasks

# target pts for predicted spline
nsamp <- length(unique(AH_long$AgeMa))
age_max <- max(AH_long$AgeMa)
age_min <- min(AH_long$AgeMa)
x_targ <- seq(age_min, age_max, by = 0.005)
head(x_targ); tail(x_targ)
nfits <- length(x_targ)
nfits

# weights
weights <- 1/(AH_long$MOT_abs_err^2)
ones = rep(1, length(AH_long$AgeMa))

# bootstrap estimates of fitted values
set.seed(1234)
nreps <- 10000

# 3 curves: Shackleton et ak, GCV optimized smoothing spline, OLS
ncurves <- 3

# arrays and vectors
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

# Reproduce the Matlab smoothing spline with the R smoothing spline

# brute-force optinmization of spar value
npts <- 201
spar_opt <- as.numeric(seq(0.0, 2.0, by = 0.01))
df_opt <- as.numeric(rep(0, npts))
lambda_opt <- as.numeric(rep(0, npts))
rss_opt <- as.numeric(rep(0, npts))
tss_opt <- as.numeric(rep(0, npts))
rsq_opt <- as.numeric(rep(0, npts))

for (k in 1:npts) {
  curve_01 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones, spar = spar_opt[k])
  df_opt[k] <- curve_01$df; lambda_opt[k] <- curve_01$lambda
  
  y_pred <- predict(curve_01, AH_long_spline_paper$AgeMa, deriv = 0)$y
  resid <- y_pred - (AH_long_spline_paper$MOT_anomaly)
  rss_opt[k] <- sum(resid^2 )
  tss_opt[k] <- sum((AH_long_spline_paper$MOT_anomaly - mean(AH_long_spline_paper$MOT_anomaly))^2)
  rsq_opt[k] <- 1.0 - (rss_opt[k]/tss_opt[k])
}

match_paper <- data.frame(spar_opt, df_opt, lambda_opt, rss_opt, tss_opt, rsq_opt)
k_opt <- which.max(rsq_opt)
print(k_opt)
print(match_paper[k_opt, ])

# Models

# Shackleton et al. spline

spar <- 0.69
curve_01 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones, spar = spar, keep.stuff = TRUE) # Shackleton et al. specs. for lambda
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

# bootsrap CL's -- calculate but do not plot
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
    # if (i < 500) {
    #   lines(x_targ, y_pred, col=adjustcolor(model_colors[j], alpha = 0.05))
    #   # points(AH_long$AgeMa, y_pred, pch = 16, cex = 0.5, col = "red")
    # }
    k <- k + 1
    fits_bootstrap[1:nfits, k] <- y_pred[1:nfits]
  }
}

boot_cl_01 <- rowQuantiles(fits_bootstrap, probs = c(0.025, 0.16, 0.50, 0.84, 0.975))
head(boot_cl_01)

# add CL's to plot
lines(x_targ, boot_cl_01[, 2], col = model_colors[1])
lines(x_targ, boot_cl_01[, 4], col = model_colors[1])
lines(x_targ, boot_cl_01[, 1], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 5], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 3], lwd = 2, col = model_colors[1])
lines(x_targ, y_pred_curve_01, lwd=2, col=model_colors[1])

#dev.print(pdf, "AH_long_Shackleton_spline.pdf")

# replot with first 500 replicates shown

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age (Ma)", ylab = "MOT anomaly (°C)", 
     font.main = 1, main = "Allen Hills (long) + 0.5 °C")
minor.tick(nx = 10, tick.ratio = 0.5)

# bootstrap estimates of fitted values
set.seed(1234)
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

lines(x_targ, boot_cl_01[, 2], col = model_colors[1])
lines(x_targ, boot_cl_01[, 4], col = model_colors[1])
lines(x_targ, boot_cl_01[, 1], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 5], lty = "dotted", lwd = 2, col = model_colors[1])
lines(x_targ, boot_cl_01[, 3], lty = "dashed", col = model_colors[1])
lines(x_targ, y_pred_curve_01, lwd=2, col=model_colors[1])

lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_anomaly, lwd=2, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_lower, lwd=1, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_upper, lwd=1, col="blue")
leg_label <- paste("Smoothing Spline, spar = ", as.character(round(spar, 3)), ", lambda = ", as.character(round(lambda[1],7)),
                   ", df = ", as.character(round(df[1], 2)), ", AIC = ", as.character(round(AIC[1], 2)),sep="")
legend(3.5, -2.75, legend=c(leg_label), col=c(model_colors[1]), lwd=2, bty="n")

points(AH_long$MOT_anomaly ~ AH_long$AgeMa, pch = 16, cex = 0.8, col = "blue")

#dev.print(pdf, "AH_long_spline_comparison.pdf")

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

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age (Ma)", ylab = "MOT anomaly (°C)", 
     font.main = 1, main = "Allen Hills (long) + 0.5 °C")
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
#dev.print(pdf, "AH_long_estimated_spline.pdf")

# OLS

curve_03 <- lm(AH_long$MOT_anomaly ~ AH_long$AgeMa) # OLS
summary(curve_03)
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

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age (Ma)", ylab = "MOT anomaly (°C)", 
     font.main = 1, main = "Allen Hills (long) + 0.5 °C")
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
