setwd("~/Dropbox/cary_projects/bean_rebuttal/bean_response/statistics/data/")
library("drc")

# activity at different pH

df <- read.csv("DODA_activity_diff_pH.csv", header = T, stringsAsFactors = T)

# DODAa1

DODAa1 <- df[ which( df$sequence == "DODAa1"), ]
DODAa1$norm_activity <- DODAa1$activity / max(DODAa1$activity) * 100
DODAa1_mean <- aggregate(activity ~ pH, data = DODAa1, FUN = mean)
DODAa1_norm_mean <- aggregate(norm_activity ~ pH, data = DODAa1, FUN = mean)

plot(activity ~ pH, data = DODAa1, col = "#FFCC00", pch = 16, xlab = "pH",
     ylab = expression(paste("Activity (", mu, "M/min)")))
lines(activity ~ pH, data = DODAa1_mean, col = "#FFCC00")
title(main = expression(paste("DODA", alpha, "1 Activity")), adj = 0)

# DODAa2

DODAa2 <- df[ which( df$sequence == "DODAa2"), ]
DODAa2$norm_activity <- DODAa2$activity / max(DODAa2$activity) * 100
DODAa2_mean <- aggregate(activity ~ pH, data = DODAa2, FUN = mean)
DODAa2_norm_mean <- aggregate(norm_activity ~ pH, data = DODAa2, FUN = mean)

plot(activity ~ pH, data = DODAa2, ylim = c(-0.1,1.3), col = "black", pch = 16, xlab = "pH",
     ylab = expression(paste("Activity (", mu, "M/min)")))
lines(activity ~ pH, data = DODAa2_mean, col = "black")
title(main = expression(paste("DODA", alpha, "2 Activity")), adj = 0)

# DODAa2-mut3

DODAa2.mut3 <- df[ which( df$sequence == "DODAa2-mut3"), ]
DODAa2.mut3$norm_activity <- DODAa2.mut3$activity / max(DODAa2.mut3$activity) * 100
DODAa2.mut3_mean <- aggregate(activity ~ pH, data = DODAa2.mut3, FUN = mean)
DODAa2.mut3_norm_mean <- aggregate(norm_activity ~ pH, data = DODAa2.mut3, FUN = mean)

plot(activity ~ pH, data = DODAa2.mut3, ylim = c(0,2.6), col = "lightgrey", pch = 16, xlab = "pH",
     ylab = expression(paste("Activity (", mu, "M/min)")))
lines(activity ~ pH, data = DODAa2.mut3_mean, col = "lightgrey")
title(main = expression(paste("DODA", alpha, "2-mut3 Activity")), adj = 0)

# all on one

plot(activity ~ pH, data = DODAa1, col = "#FFCC00", ylim = c(0, 2.6), xlim = c(4.9, 8.6), pch = 16, xlab = "pH",
     ylab = expression(paste("Activity (", mu, "M/min)")), xaxt = "n")
axis(1, at = c(5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5), las=2)
lines(activity ~ unique(pH), data = DODAa1_mean, col = "#FFCC00")
#lines(spline(jitter_pH, DODAa1$activity, n=100), col = "#FFCC00")
points(activity ~ pH, data = DODAa2, col = "black", pch = 17)
lines(activity ~ pH, data = DODAa2_mean, col = "black")
#lines(spline(DODAa2$pH, DODAa2$activity, n=100), col = "darkgrey")
points(activity ~ pH, data = DODAa2.mut3, col = "lightgrey", pch = 18)
lines(activity ~ unique(pH), data = DODAa2.mut3_mean, col = "lightgrey")
#lines(spline(jitter_pH, DODAa2.mut3$activity, n = 100), col = "lightgrey")
legend("topleft", legend = c(expression(paste("DODA", alpha, "1")), expression(paste("DODA", alpha, "2")), 
                             expression(paste("DODA", alpha, "2-mut3"))),
       pch = 16:18, col = c("#FFCC00", "black", "lightgrey"))
title(main = expression(paste("DODA activity at differing pH")), adj=0)

# all on one normalised

plot(norm_activity ~ pH, data = DODAa1, col = "#FFCC00", ylim = c(0, 100), xlim = c(4.9, 8.6), pch = 16, xlab = "pH",
     ylab = expression(paste("Normalised activity (%)")), xaxt = "n")
axis(1, at = c(5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5), las=2)
lines(norm_activity ~ pH, data = DODAa1_norm_mean, col = "#FFCC00")
points(norm_activity ~ pH, data = DODAa2, col = "black", pch = 17)
lines(norm_activity ~ pH, data = DODAa2_norm_mean, col = "black")
points(norm_activity ~ pH, data = DODAa2.mut3, col = "lightgrey", pch = 18)
lines(norm_activity ~ pH, data = DODAa2.mut3_norm_mean, col = "lightgrey")
legend("topleft", legend = c(expression(paste("DODA", alpha, "1")), expression(paste("DODA", alpha, "2")), 
                             expression(paste("DODA", alpha, "2-mut3"))),
       pch = 16:18, col = c("#FFCC00", "black", "lightgrey"))
title(main = expression(paste("Normalised DODA activity at differing pH")), adj=0)

# pH 8.5
df <- read.csv("DODAa2_a2-mut3_activity_pH8.5.csv", header = T, stringsAsFactors = T)

# DODAa2

DODAa2 <- df[ which( df$sequence == "DODAa2"), ]
DODAa2_mean <- aggregate(activity ~ ldopa, data = DODAa2, FUN = mean)

start <- c(5,2)
names(start) <- c("a", "b")

fit.mm <- nls(activity ~ a * ldopa / (b + ldopa), data = DODAa2, start = start)
fit.drm <- drm(activity ~ ldopa, data = DODAa2, fct = MM.2())
summary(fit.mm)

# Formula: activity ~ a * ldopa/(b + ldopa)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a   4.3049     0.5263   8.179 8.27e-08 ***
# b   1.8935     0.4932   3.839  0.00102 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2797 on 20 degrees of freedom
# 
# Number of iterations to convergence: 8 
# Achieved convergence tolerance: 4.053e-06

summary(fit.drm)

# Model fitted: Michaelis-Menten (2 parms)
# 
# Parameter estimates:
#   
#   Estimate Std. Error t-value   p-value    
# d:(Intercept)  4.30486    0.62994  6.8338 1.214e-06 ***
# e:(Intercept)  1.89350    0.59818  3.1654  0.004865 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error:
#   
#   0.2797152 (20 degrees of freedom)

mm <- function(ldopa, a, b) { a * ldopa / (b + ldopa) }

plot(activity ~ ldopa, data=DODAa2, ylim = c(0,3.5), col = "black", pch = 16, ylab = expression(paste("Activity (", mu, "M/min)")),
     xlab = expression(paste("[L-DOPA (mM)]")))
lines(seq(0,3.8,0.01), mm(seq(0,3.8,0.01), 4.30486, 1.8935), col = "black") # using estimated params
title(main = expression(paste("DODA", alpha, "2, pH 8.5")), adj = 0)
text(3.5, 1.0, labels = expression(paste("V"["M"] * " = 4.305")))
text(3.5, 0.8, labels = expression(paste("K"["M"] * " = 1.894")))

# DODAa2-mut3

DODAa2.mut3 <- df[ which( df$sequence == "DODAa2-mut3"), ]

start <- c(10,4)
names(start) <- c("a", "b")

fit.mm <- nls(activity ~ a * ldopa / (b + ldopa), data = DODAa2.mut3, start = start)
fit.drm <- drm(activity ~ ldopa, data = DODAa2.mut3, fct = MM.2())

summary(fit.mm)

# Formula: activity ~ a * ldopa/(b + ldopa)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a   6.5110     0.4182  15.567 1.21e-12 ***
# b   1.2279     0.1970   6.234 4.35e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3463 on 20 degrees of freedom
# 
# Number of iterations to convergence: 7 
# Achieved convergence tolerance: 5.449e-06

summary(fit.drm)

# Model fitted: Michaelis-Menten (2 parms)
# 
# Parameter estimates:
#   
#   Estimate Std. Error t-value   p-value    
# d:(Intercept)  6.51099    0.42860 15.1913 1.905e-12 ***
# e:(Intercept)  1.22791    0.20261  6.0605 6.344e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error:
#   
#   0.3463209 (20 degrees of freedom)

plot(activity ~ ldopa, data=DODAa2.mut3, ylim = c(0,5.5), col = "lightgrey", pch = 16, ylab = expression(paste("Activity (", mu, "M/min)")),
     xlab = expression(paste("[L-DOPA (mM)]")))
lines(seq(0,3.8,0.01), mm(seq(0,3.8,0.01), 6.51099, 1.22791), col = "lightgrey") # using estimated params
title(main = expression(paste("DODA", alpha, "2-mut3, pH 8.5")), adj = 0)
text(3.5, 2.0, labels = expression(paste("V"["M"] * " = 6.511")))
text(3.5, 1.8, labels = expression(paste("K"["M"] * " = 1.228")))

# DODAa2-mut3 and DODAa2 plotted together

plot(activity ~ ldopa, data=DODAa2.mut3, ylim = c(0,5.5), col = "lightgrey", pch = 16, ylab = expression(paste("Activity (", mu, "M/min)")),
     xlab = expression(paste("[L-DOPA (mM)]")))
lines(seq(0,3.8,0.01), mm(seq(0,3.8,0.01), 6.51099, 1.22791), col = "lightgrey")

text(0.5, 5.0, labels = expression(paste("V"["M"] * " = 6.511")))
text(0.5, 4.7, labels = expression(paste("K"["M"] * " = 1.228")))

points(DODAa2$ldopa, DODAa2$activity, col = "black", pch=17)
lines(seq(0,3.8,0.01), mm(seq(0,3.8,0.01), 4.30486, 1.8935), col = "black")

text(3.5, 1.5, labels = expression(paste("V"["M"] * " = 4.305")))
text(3.5, 1.2, labels = expression(paste("K"["M"] * " = 1.894")))

title(main = expression(paste("DODA", alpha, "2-mut3 and DODA", alpha, "2, pH 8.5")), adj = 0)
legend("bottom", legend = c(expression(paste("DODA", alpha, "2-mut3")), expression(paste("DODA", alpha, "2"))),
       pch = 16:17, col = c("lightgrey", "black"))

# pH 6.0

df <- read.csv("DODA_activity_pH6.csv", header = T, stringsAsFactors = T)

# DODAa1

DODAa1 <- df[ which( df$sequence == "DODAa1"), ]

start <- c(300, 0.5, 1)
names(start) <- c("a", "b", "c")

fit.inhib <- nls(activity ~ a * ldopa / (b + ldopa + ldopa^2 / c), data = DODAa1, start = start)
summary(fit.inhib)

# Formula: activity ~ a * ldopa/(b + ldopa + ldopa^2/c)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)
# a  47.65820  158.37197   0.301    0.767
# b   2.73047    9.55832   0.286    0.778
# c   0.02808    0.09656   0.291    0.774
# 
# Residual standard error: 0.1851 on 19 degrees of freedom
# 
# Number of iterations to convergence: 24 
# Achieved convergence tolerance: 4.859e-06

inhib <- function(ldopa, a, b, c) { a * ldopa / (b + ldopa + ldopa^2 / c) }

plot(activity ~ ldopa, data=DODAa1, ylim = c(0,2.6), col = "#FFCC00", pch = 16, ylab = expression(paste("Activity (", mu, "M/min)")),
     xlab = expression(paste("[L-DOPA (mM)]")))
lines(seq(0,3.8,0.01), inhib(seq(0,3.8,0.01), 47.65820, 2.73047, 0.02808), col = "#FFCC00") # using estimated params
title(main = expression(paste("DODA", alpha, "1, pH 6.0")), adj = 0)
text(3.0, 2.5, labels = expression(paste("V"["M"] * " = 47.658")))
text(3.0, 2.35, labels = expression(paste("K"["M"] * " = 2.73")))
text(3.0, 2.2, labels = expression(paste("K"["I"] * " = 0.0281")))

# DODAa2

DODAa2 <- df[ which( df$sequence == "DODAa2"), ]

start <- c(5,2)
names(start) <- c("a", "b")

fit.mm <- nls(activity ~ a * ldopa / (b + ldopa), data = DODAa2, start = start)
summary(fit.mm)
# 
# Formula: activity ~ a * ldopa/(b + ldopa)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a  0.44639    0.02074  21.520 2.66e-15 ***
# b  0.03009    0.02013   1.495    0.151    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06367 on 20 degrees of freedom
# 
# Number of iterations to convergence: 12 
# Achieved convergence tolerance: 3.325e-06

plot(activity ~ ldopa, data=DODAa2, ylim = c(0,0.6), col = "black", pch = 16, ylab = expression(paste("Activity (", mu, "M/min)")),
     xlab = expression(paste("[L-DOPA (mM)]")))
lines(seq(0,3.8,0.01), mm(seq(0,3.8,0.01), 0.44639, 0.03009), col = "black") # using estimated params
title(main = expression(paste("DODA", alpha, "2, pH 6.0")), adj = 0)
text(3.5, 0.2, labels = expression(paste("V"["M"] * " = 0.446")))
text(3.5, 0.17, labels = expression(paste("K"["M"] * " = 0.03")))

# DODAa2-mut3

DODAa2.mut3 <- df[ which( df$sequence == "DODAa2-mut3"), ]

start <- c(5,2)
names(start) <- c("a", "b")

fit.mm <- nls(activity ~ a * ldopa / (b + ldopa), data = DODAa2.mut3, start = start)
summary(fit.mm)

# Formula: activity ~ a * ldopa/(b + ldopa)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a  0.48195    0.02972  16.214 5.68e-13 ***
# b  0.08184    0.03719   2.201   0.0397 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08133 on 20 degrees of freedom
# 
# Number of iterations to convergence: 12 
# Achieved convergence tolerance: 4.281e-06

plot(activity ~ ldopa, data=DODAa2.mut3, ylim = c(0,0.62), col = "lightgrey", pch = 16, ylab = expression(paste("Activity (", mu, "M/min)")),
     xlab = expression(paste("[L-DOPA (mM)]")))
lines(seq(0,3.8,0.01), mm(seq(0,3.8,0.01), 0.48195, 0.08184), col = "lightgrey") # using estimated params
title(main = expression(paste("DODA", alpha, "2-mut3, pH 6.0")), adj = 0)
text(3.5, 0.2, labels = expression(paste("V"["M"] * " = 0.482")))
text(3.5, 0.17, labels = expression(paste("K"["M"] * " = 0.082")))

# all three on one plot at pH 6.0

plot(activity ~ ldopa, data=DODAa2.mut3, ylim = c(0,2.6), col = "lightgrey", pch = 18, ylab = expression(paste("Activity (", mu, "M/min)")),
     xlab = expression(paste("[L-DOPA (mM)]")))
points(DODAa2$ldopa, DODAa2$activity, col = "black", pch=17)
points(DODAa1$ldopa, DODAa1$activity, col = "#FFCC00", pch=16)
lines(seq(0,3.8,0.01), inhib(seq(0,3.8,0.01), 47.65820, 2.73047, 0.02808), col = "#FFCC00")
lines(seq(0,3.8,0.01), mm(seq(0,3.8,0.01), 0.44639, 0.03009), col = "black")
lines(seq(0,3.8,0.01), mm(seq(0,3.8,0.01), 0.48195, 0.08184), col = "lightgrey")
legend("topright", legend = c(expression(paste("DODA", alpha, "1")), expression(paste("DODA", alpha, "2")), 
                             expression(paste("DODA", alpha, "2-mut3"))),
       pch = 16:18, col = c("#FFCC00", "black", "lightgrey"))
title(main = expression(paste("DODA kinetics at pH 6.0")), adj = 0)
