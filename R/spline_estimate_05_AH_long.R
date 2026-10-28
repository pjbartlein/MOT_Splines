library(readxl, writexl)
library(Hmisc)
library(matrixStats)

# bootstraped splines =============================================================================================

# Allen_Hills_long
data_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/MOTs/Shackleton/"
file_name <- "AH_long.xlsx"
xlxs_file <- paste(data_path, file_name, sep = "")
AH_long <- read_excel(xlxs_file)
head(AH_long); tail(AH_long)

# MOT + 0.5C
AH_long$MOT_anomaly <- AH_long$MOT_anomaly + 0.5

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

# set up

weights <- 1/(AH_long$MOT_abs_err^2)
ones = rep(1, length(AH_long$AgeMa))

# bootstrap estimates of fitted values
set.seed(1234)

ncurves <- 3
nreps <- 10000

fits <- matrix(0, nrow = nsamp, ncol = ncurves)
residuals <- matrix(0, nrow = nsamp, ncol = ncurves)
fits_bootstrap <- matrix(0, nrow=nsamp, ncol = nreps)
df <- as.numeric(rep(0, ncurves))
lambda <- as.numeric(rep(0, ncurves))
AIC <- as.numeric(rep(0, ncurves))
rsq <- as.numeric(rep(0, ncurves))
adj_rsq <- as.numeric(rep(0, ncurves))

model_colors <- c("black", "purple", "red")

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT_anomaly", main = "Allen_Hills_long + 0.5 C")

# add spline from paper
file_name <- "AH_long_spline_paper.xlsx"
xlxs_file <- paste(data_path, file_name, sep = "")
AH_long_spline_paper <- read_excel(xlxs_file)
head(AH_long_spline_paper); tail(AH_long_spline_paper)

lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_anomaly + 0.5, lwd=2, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_lower + 0.5, lwd=1, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_upper + 0.5, lwd=1, col="blue")
leg_label <- paste("Spline from Shackleton et al. 2026 (Fig. 1)", sep="")
legend(3.5, -2.75, legend=c(leg_label), col="blue", lwd=2, bty="n")
# dev.print(pdf, "AH_long_Shackleton_spline.pdf")

# Models

# Shackleton et al. spline

spar <- 0.0
for (k in 1:201) {
  # curve_01 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones, lambda = lambda1) # Shackleton et al. specs. for lambda
  # curve_01 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones, df = 13) # Best guess
  curve_01 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones, spar = spar) # Best guess
  df[1] <- curve_01$df; lambda[1] <- curve_01$lambda; fits[, 1] <- curve_01$y
  residuals[, 1] <- AH_long$MOT_anomaly - curve_01$y
  sigma2 <- sum(residuals[, 1]^2)/(nsamp - 1.0)
  rss <- sum(residuals[, 1]^2 )
  tss <- sum((AH_long$MOT_anomaly - mean(AH_long$MOT_anomaly))^2)
  # http://users.stat.umn.edu/~helwig/notes/smooth-spline-notes.html#aic-and-bic
  max_log_lik <- (-1.0 / (2.0 * sigma2)) * rss - 1.0 * (nsamp / 2.0) * log(sigma2) - 1.0 * (nsamp / 2.0) * log(2 * pi)
  AIC[1] <- - 2.0 * max_log_lik + 2.0 * df[1]
  y_pred <- predict(curve_01, AH_long_spline_paper$AgeMa, deriv = 0)$y
  resid <- y_pred - (AH_long_spline_paper$MOT_anomaly + 0.5)
  rss_spline <- sum(resid^2 )
  tss_spline <- sum((AH_long_spline_paper$MOT_anomaly + 0.5) - mean(AH_long_spline_paper$MOT_anomaly) + 0.5)
  rsq_spline <- 1.0 - (rss_spline/tss_spline)
  # AIC[1] <- round(AIC[1], 2)
  # print(AIC[1])
  
  print (c(spar, df[1], lambda[1], rss, AIC[1], rss_spline, tss_spline, rsq_spline))
  spar <- spar + 0.01
}
  
  # R-sq values
  rsq[1] <- 1 - (rss/tss)
  adj_rsq[1] <- 1 - (rss / (nsamp - df[1])) / (tss / (nsamp - 1))
  rsq[1] <- round(rsq[1], digits = 3)
  adj_rsq[1] <- round(adj_rsq[1], digits = 3)
  rsq[1]; adj_rsq[1]

hist(residuals[, 1], breaks = seq(-4.25, 4.25, by = 0.25))
qqnorm(residuals[, 1]); qqline(residuals[, 1])
shapiro.test(residuals[, 1])

oldpar <- par(pty="s")
ylim_resid <- c(-2, 2)
plot(residuals[, 1] ~ fits[, 1], ylim = ylim_resid, xlim = ylim_resid, pch = 16, cex = 0.8, 
     ylab = "Residuals from Smoothing Spline", xlab =  "Smoothing Spline Fitted Values")
lo_resid <- loess(residuals[, 1] ~ fits[, 1], span=0.80, degree = 1, iterations = 3)
lines(lo_resid$fitted ~ fits[, 1], col = "red", lwd = 2)
par(oldpar)

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT_anomaly", main = "Allen_Hills_long + 0.5 C")

# bootsrap CL's
fits_bootstrap <- matrix(0, nrow=nsamp, ncol = nreps)
j <- 1
k <- 0
for (i in 1:nreps) {
  # print(i)
  x_boot <- AH_long$AgeMa + rnorm(nsamp, mean = 0, sd = AH_long$Rel_err_Myr)
  x_boot <- x_boot[order(x_boot)]
  y_boot <- AH_long$MOT_anomaly + rnorm(nsamp, mean = 0, sd = AH_long$MOT_abs_err)
  y_boot <- y_boot[order(x_boot)]
  # points(x_boot, y_boot, pch = 16, cex = 0.5, col = "red")
  
  fit_boot <- smooth.spline(x_boot, y_boot, w = ones, df = df[1])
  if (length(x_boot) == length(fit_boot$y)) {
    y_pred <- predict(fit_boot, AH_long$AgeMa, deriv = 0)$y
    if (i < 1000) {
      lines(AH_long$AgeMa, y_pred, col=adjustcolor(model_colors[j], alpha = 0.01))
      # points(AH_long$AgeMa, y_pred, pch = 16, cex = 0.5, col = "red")
    }
    k <- k + 1
    fits_bootstrap[1:nsamp, k] <- y_pred[1:nsamp]
  }
}

boot_cl_01 <- rowQuantiles(fits_bootstrap, probs = c(0.025, 0.16, 0.50, 0.84, 0.975))
head(boot_cl_01)
lines(AH_long$AgeMa, boot_cl_01[, 2], col = model_colors[1])
lines(AH_long$AgeMa, boot_cl_01[, 4], col = model_colors[1])
lines(AH_long$AgeMa, boot_cl_01[, 1], lty = "dotted", col = model_colors[1])
lines(AH_long$AgeMa, boot_cl_01[, 5], lty = "dotted", col = model_colors[1])
lines(AH_long$AgeMa, boot_cl_01[, 3], lty = "dashed", col = model_colors[1])

# Shackleton et al. specs. for lambda
lines(AH_long$AgeMa, fits[, 1], lwd=2, col=model_colors[1])
leg_label <- paste("Smoothing Spline, lambda = ", as.character(round(lambda[1],7)), 
                   ", df = ", as.character(round(df[1], 2)), ", AIC = ", as.character(round(AIC[1], 2)),sep="")
legend(3.5, -2.75, legend=c(leg_label), col=c(model_colors[1]), lwd=2, bty="n")

# spline from paper
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_anomaly + 0.5, lwd=2, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_lower + 0.5, lwd=1, col="blue")
lines(AH_long_spline_paper$AgeMa, AH_long_spline_paper$MOT_upper + 0.5, lwd=1, col="blue")
leg_label <- paste("Spline from Shackleton et al. 2026 (Fig. 1)", sep="")
legend(3.5, -2.50, legend=c(leg_label), col="blue", lwd=2, bty="n")

points(AH_long$MOT_anomaly ~ AH_long$AgeMa, pch = 16, cex = 0.8, col = "blue")

dev.print(pdf, "AH_long_Shackleton_spline.pdf")


# optimal

curve_02 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = ones) # estimate lambda, spar, df
df[2] <- curve_02$df; lambda[2] <- curve_02$lambda; fits[, 2] <- curve_02$y
residuals[, 2] <- AH_long$MOT_anomaly - curve_02$y
sigma2 <- sum(residuals[, 2]^2)/(nsamp - 1.0)
rss <- sum(residuals[, 2]^2 )
tss <- sum((AH_long$MOT_anomaly - mean(AH_long$MOT_anomaly))^2)
# http://users.stat.umn.edu/~helwig/notes/smooth-spline-notes.html#aic-and-bic
max_log_lik <- (-1.0 / (2.0 * sigma2)) * rss - 1.0 * (nsamp / 2.0) * log(sigma2) - 1.0 * (nsamp / 2.0) * log(2 * pi)
AIC[2] <- - 2.0 * max_log_lik + 2.0 * df[2]
AIC[2] <- round(AIC[2], 2)
print(AIC[2])

# R-sq values
rsq[2] <- 1 - (rss/tss)
adj_rsq[2] <- 1 - (rss / (nsamp - df[2])) / (tss / (nsamp - 1))
rsq[2] <- round(rsq[2], digits = 3)
adj_rsq[2] <- round(adj_rsq[2], digits = 3)
rsq[2]; adj_rsq[2]

hist(residuals[, 2], breaks = seq(-4.25, 4.25, by = 0.25))
qqnorm(residuals[, 2]); qqline(residuals[, 2])
shapiro.test(residuals[, 2])

oldpar <- par(pty="s")
ylim_resid <- c(-2, 2)
plot(residuals[, 2] ~ fits[, 2], ylim = ylim_resid, xlim = ylim_resid, pch = 16, cex = 0.8, 
    ylab = "Residuals from Smoothing Spline", xlab =  "Smoothing Spline Fitted Values")
lo_resid <- loess(residuals[, 2] ~ fits[, 2], span=0.80, degree = 1, iterations = 3)
lines(lo_resid$fitted ~ fits[, 2], col = "red", lwd = 2)
par(oldpar)

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT_anomaly", main = "Allen_Hills_long + 0.5 C")

# bootsrap CL's
fits_bootstrap <- matrix(0, nrow=nsamp, ncol = nreps)
k <- 0
for (i in 1:nreps) {
  # print(i)
  x_boot <- AH_long$AgeMa + rnorm(nsamp, mean = 0, sd = AH_long$Rel_err_Myr)
  x_boot <- x_boot[order(x_boot)]
  y_boot <- AH_long$MOT_anomaly + rnorm(nsamp, mean = 0, sd = AH_long$MOT_abs_err)
  y_boot <- y_boot[order(x_boot)]
  # points(x_boot, y_boot, pch = 16, cex = 0.5, col = "red")
  
  fit_boot <- smooth.spline(x_boot, y_boot, w = ones, df = df[2])
  if (length(x_boot) == length(fit_boot$y)) {
    y_pred <- predict(fit_boot, AH_long$AgeMa, deriv = 0)$y
    if (i < 1000) {
      lines(AH_long$AgeMa, y_pred, col=adjustcolor(model_colors[2], alpha = 0.01))
      # points(AH_long$AgeMa, y_pred, pch = 16, cex = 0.5, col = "red")
    }
    k <- k + 1
    fits_bootstrap[1:nsamp, k] <- y_pred[1:nsamp]
  }
}

boot_cl_02 <- rowQuantiles(fits_bootstrap, probs = c(0.025, 0.16, 0.50, 0.84, 0.975))
head(boot_cl_01)
lines(AH_long$AgeMa, boot_cl_02[, 2], col = model_colors[2])
lines(AH_long$AgeMa, boot_cl_02[, 4], col = model_colors[2])
lines(AH_long$AgeMa, boot_cl_02[, 1], lty = "dotted", col = model_colors[2])
lines(AH_long$AgeMa, boot_cl_02[, 5], lty = "dotted", col = model_colors[2])
lines(AH_long$AgeMa, boot_cl_02[, 3], lty = "dashed", col = model_colors[2])

# optimal lambda
lines(AH_long$AgeMa, fits[, 2], lwd=2, col=model_colors[2])
leg_label <- paste("Smoothing Spline, estimated lambda = ", as.character(round(lambda[2],7)), 
                   ", df = ", as.character(round(df[2], 2)), ", AIC = ", as.character(round(AIC[2], 2)),sep="")
legend(3.5, -2.75, legend=c(leg_label), col=c(model_colors[2]), lwd=2, bty="n")

points(AH_long$MOT_anomaly ~ AH_long$AgeMa, pch = 16, cex = 0.8, col = "blue")
dev.print(pdf, "AH_long_estimated_spline.pdf")

# OLS

curve_03 <- lm(AH_long$MOT_anomaly ~ AH_long$AgeMa) # OLS
fits[, 3] <- curve_03$fitted.values
residuals[, 3] <- AH_long$MOT_anomaly - curve_03$fitted.values
rss <- sum(residuals[, 3]^2 )
tss <- sum((AH_long$MOT_anomaly - mean(AH_long$MOT_anomaly))^2)
AIC[3] <- round(AIC(curve_03), 2)
print(AIC[3])

# R-sq values
rsq[3] <- 1 - (rss/tss)
adj_rsq[3] <- 1 - (rss / (nsamp - df[3])) / (tss / (nsamp - 1))
rsq[3] <- round(rsq[3], digits = 3)
adj_rsq[3] <- round(adj_rsq[3], digits = 3)
rsq[3]; adj_rsq[3]

hist(residuals[, 3], breaks = seq(-4.25, 4.25, by = 0.25))
qqnorm(residuals[, 3]); qqline(residuals[, 3])
shapiro.test(residuals[, 3])

oldpar <- par(pty="s")
ylim_resid <- c(-2, 2)
plot(residuals[, 3] ~ fits[, 3], ylim = ylim_resid, xlim = ylim_resid, pch = 16, cex = 0.8, 
     ylab = "Residuals from Smoothing Spline", xlab =  "Smoothing Spline Fitted Values")
lo_resid <- loess(residuals[, 3] ~ fits[, 3], span=0.80, degree = 1, iterations = 3)
lines(lo_resid$fitted ~ fits[, 3], col = "red", lwd = 2)
par(oldpar)

plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT_anomaly", main = "Allen_Hills_long + 0.5 C")

# bootsrap CL's
fits_bootstrap <- matrix(0, nrow=nsamp, ncol = nreps)
k <- 0
for (i in 1:nreps) {
  # print(i)
  x_boot <- AH_long$AgeMa + rnorm(nsamp, mean = 0, sd = AH_long$Rel_err_Myr)
  x_boot <- x_boot[order(x_boot)]
  y_boot <- AH_long$MOT_anomaly + rnorm(nsamp, mean = 0, sd = AH_long$MOT_abs_err)
  y_boot <- y_boot[order(x_boot)]
  # points(x_boot, y_boot, pch = 16, cex = 0.5, col = "red")
  
  fit_boot <- lm(y_boot ~ x_boot) # OLS
  if (length(x_boot) == length(fit_boot$fitted.values)) {
    y_pred <- fit_boot$coefficients[1] + fit_boot$coefficients[2] * AH_long$AgeMa
    if (i < 1000) {
      lines(AH_long$AgeMa, y_pred, col=adjustcolor(model_colors[3], alpha = 0.01))
      # points(AH_long$AgeMa, y_pred, pch = 16, cex = 0.5, col = "red")
    }
    k <- k + 1
    fits_bootstrap[1:nsamp, k] <- y_pred[1:nsamp]
  }
}
boot_cl_03 <- rowQuantiles(fits_bootstrap, probs = c(0.025, 0.16, 0.50, 0.84, 0.975))
head(boot_cl_03)
lines(AH_long$AgeMa, boot_cl_03[, 2], col = model_colors[3])
lines(AH_long$AgeMa, boot_cl_03[, 4], col = model_colors[3])
lines(AH_long$AgeMa, boot_cl_03[, 1], lty = "dotted", col = model_colors[3])
lines(AH_long$AgeMa, boot_cl_03[, 5], lty = "dotted", col = model_colors[3])
lines(AH_long$AgeMa, boot_cl_03[, 3], lty = "dashed", col = model_colors[3])

lines(AH_long$AgeMa, fits[, 3], lwd=2, col=model_colors[3])
leg_label <- paste("OLS Linear Model", ", AIC = ", as.character(round(AIC[3], 2)), sep="")
legend(3.5, -2.75, legend=c(leg_label), col=model_colors[3], lwd=2, bty="n")

points(AH_long$MOT_anomaly ~ AH_long$AgeMa, pch = 16, cex = 0.8, col = "blue")
dev.print(pdf, "AH_long_OLS.pdf")


plot(AH_long$MOT_anomaly ~ AH_long$AgeMa, type = "p", pch = 16, cex = 0.8, col = "blue",
     xlim = c(3.5, 0), ylim = c(-3.0, 3.0), xlab = "Age Ma", ylab = "MOT_anomaly", main = "Allen_Hills_long + 0.5 C")

lines(AH_long$AgeMa, boot_cl_01[, 2], col = model_colors[1])
lines(AH_long$AgeMa, boot_cl_01[, 4], col = model_colors[1])
lines(AH_long$AgeMa, boot_cl_01[, 1], lty = "dotted", col = model_colors[1])
lines(AH_long$AgeMa, boot_cl_01[, 5], lty = "dotted", col = model_colors[1])
lines(AH_long$AgeMa, fits[, 1], lwd=2, col=model_colors[1])
leg_label <- paste("Smoothing Spline, lambda = ", as.character(round(lambda[1],7)), 
                   ", df = ", as.character(round(df[1], 2)), ", AIC = ", as.character(round(AIC[1], 2)),sep="")
legend(3.5, -2.25, legend=c(leg_label), col=c(model_colors[1]), lwd=2, bty="n")

lines(AH_long$AgeMa, boot_cl_02[, 2], col = model_colors[2])
lines(AH_long$AgeMa, boot_cl_02[, 4], col = model_colors[2])
lines(AH_long$AgeMa, boot_cl_02[, 1], lty = "dotted", col = model_colors[2])
lines(AH_long$AgeMa, boot_cl_02[, 5], lty = "dotted", col = model_colors[2])
lines(AH_long$AgeMa, fits[, 2], lwd=2, col=model_colors[2])
lines(AH_long$AgeMa, fits[, 2], lwd=2, col=model_colors[2])
leg_label <- paste("Smoothing Spline, estimated lambda = ", as.character(round(lambda[2],7)), 
                   ", df = ", as.character(round(df[2], 2)), ", AIC = ", as.character(round(AIC[2], 2)),sep="")
legend(3.5, -2.50, legend=c(leg_label), col=c(model_colors[2]), lwd=2, bty="n")

lines(AH_long$AgeMa, boot_cl_03[, 2], col = model_colors[3])
lines(AH_long$AgeMa, boot_cl_03[, 4], col = model_colors[3])
lines(AH_long$AgeMa, boot_cl_03[, 1], lty = "dotted", col = model_colors[3])
lines(AH_long$AgeMa, boot_cl_03[, 5], lty = "dotted", col = model_colors[3])
lines(AH_long$AgeMa, fits[, 3], lwd=2, col=model_colors[3])
leg_label <- paste("OLS Linear Model", ", AIC = ", as.character(round(AIC[3], 2)), sep="")
legend(3.5, -2.75, legend=c(leg_label), col=model_colors[3], lwd=2, bty="n")
dev.print(pdf, "AH_long_all.pdf")

out_data <- data.frame(AH_long$AgeMa, AH_long$MOT_anomaly, 
                       fits[, 1], boot_cl_01[, 1], boot_cl_01[, 2], boot_cl_01[, 3], boot_cl_01[, 4], boot_cl_01[, 5],
                       fits[, 2], boot_cl_02[, 1], boot_cl_02[, 2], boot_cl_02[, 3], boot_cl_02[, 4], boot_cl_02[, 5],
                       fits[, 3], boot_cl_03[, 1], boot_cl_03[, 2], boot_cl_03[, 3], boot_cl_03[, 4], boot_cl_03[, 5] )
head(out_data)
names(out_data) <- c("AgeMA", "MOT_anomaly", 
                     "Shack_fit", "Shack_q0.025", "Shack_q0.16", "Shack_median", "Shack_q0.84", "Shack_q0.975",
                     "opt_fit", "opt_q0.025", "opt_q0.16", "opt_median", "opt_q0.84", "opt_q0.975",
                     "linear_fit", "linear_q0.025", "linear_q0.16", "linear_median", "linear_q0.84", "linear_q0.975" )
head(out_data)
write.csv(out_data, "AH+0.5_splines.csv", row.names = FALSE)

# # ad-hoc approx
# df = T / P
# spline_AH_long_02 <- smooth.spline(AH_long$AgeMa, AH_long$MOT_anomaly, w = weights, df = df)
# spline_AH_long_02
# lines(spline_AH_long_02$x, spline_AH_long_02$y, lwd=2, col="lightblue")
# leg_label <- paste("Weighted Smoothing Spline, df = T/P = ", as.character(round(spline_AH_long_02$df,2)), 
#                    ", lambda = ", as.character(round(spline_AH_long_02$lambda,4)), sep="")
# legend(3.5, -3.25, legend=c(leg_label), col=c("lightblue"), lwd=2, bty="n")


