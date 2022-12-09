require(ggplot2)

setwd("~/Dropbox/cary_projects/bean_rebuttal/bean_response/statistics/data/")
#setwd("D:/Dropbox/cary_projects/DODA/cary/bean_rebuttal/bean_response/statistics/data/")

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# Nicotiana heterologous, our (MN) BvDODAa2
df1 <- read.csv("Nbenthamiana_MN.csv", header=T, stringsAsFactors = T)
# this is redone with no background subtraction, so slightly different to previous values
str(df1) 
df1$sequence <- factor(df1$sequence, c("DODAa1", "DODAa2-MN-mut3", "DODAa2-MN", "uninfiltrated"))

# plotting

means <- aggregate(betanin_ngmg ~ sequence, data = df1, FUN = mean)
rownames(means) <- means[,1]
means <- means[,-1]
sds <- aggregate(betanin_ngmg ~ sequence, data = df1, FUN = sd)
rownames(sds) <- sds[,1]
sds <- sds[,-1]

bp <- barplot(means, names.arg = c(expression(paste("DODA", alpha, "1")),
                                   expression(paste("DODA", alpha, "2-mut3")), 
                                   expression(paste("DODA", alpha, "2")),
                                   expression(paste("Control"))),
              col = "#a41389", border = "NA", yaxt = "n", ylim = c(-5, 300))

error.bar(bp, means, sds)

axis(2, at = c(0, 50, 100, 150, 200, 250, 300))
title(main = expression(paste("DODA expression in ", italic("N. benthamiana"), " leaves")), adj = 0)
title(ylab = expression(paste("Betanin (ng/mg FW)")))

points(jitter(rep(.5, 5)), df1$betanin_ngmg[ which(df1$sequence == "DODAa1") ], pch = 1)
points(jitter(rep(1.7, 5)), df1$betanin_ngmg[ which(df1$sequence == "DODAa2-MN-mut3") ], pch = 1)
points(jitter(rep(2.9, 3)), df1$betanin_ngmg[ which(df1$sequence == "DODAa2-MN") ], pch = 1)
points(jitter(rep(4.1, 4)), df1$betanin_ngmg[ which(df1$sequence == "uninfiltrated") ], pch = 1)


X <- grconvertX(c(1.5, 5), from = "user", to = "ndc")
Y <- grconvertY(c(25, 300), from = "user", to = "ndc")

par(fig = c(X, Y), new = T)

bp <- barplot(means[-1], col = "#a41389", border = "NA", yaxt = "n", ylim = c(-0.1, 5))
error.bar(bp, means[-1], sds[-1])
axis(2, at = c(0, 1, 2, 3), labels = c(0, 1, 2, 3))
points(jitter(rep(.5, 5)), df1$betanin_ngmg[ which(df1$sequence == "DODAa2-MN-mut3") ], pch = 1)
points(jitter(rep(1.7, 3)), df1$betanin_ngmg[ which(df1$sequence == "DODAa2-MN") ], pch = 1)
points(jitter(rep(2.9, 4)), df1$betanin_ngmg[ which(df1$sequence == "uninfiltrated") ], pch = 1)

# stats

t.test(df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN-mut3")], df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN")],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN-mut3")] and df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN")]
# t = -1.4984, df = 2.5164, p-value = 0.8762
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   -2.703253       Inf
# sample estimates:
#   mean of x mean of y 
# 0.9363474 1.9369473 

t.test(df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN-mut3")], df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN")],
       var.equal = F, alternative = "two.sided")

# Welch Two Sample t-test
# 
# data:  df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN-mut3")] and df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN")]
# t = -1.4984, df = 2.5164, p-value = 0.2476
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -3.376608  1.375409
# sample estimates:
#   mean of x mean of y 
# 0.9363474 1.9369473 

t.test(df1$betanin_ngmg[which(df1$sequence == "DODAa1")], df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN")],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df1$betanin_ngmg[which(df1$sequence == "DODAa1")] and df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN")]
# t = 6.2007, df = 4.0026, p-value = 0.001717
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   142.6554      Inf
# sample estimates:
#   mean of x  mean of y 
# 219.314675   1.936947 

t.test(df1$betanin_ngmg[which(df1$sequence == "DODAa1")], df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN-mut3")],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df1$betanin_ngmg[which(df1$sequence == "DODAa1")] and df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN-mut3")]
# t = 6.2301, df = 4.0003, p-value = 0.00169
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   143.6545      Inf
# sample estimates:
#   mean of x   mean of y 
# 219.3146750   0.9363474 

t.test(df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN")], df1$betanin_ngmg[which(df1$sequence == "uninfiltrated")],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN")] and df1$betanin_ngmg[which(df1$sequence == "uninfiltrated")]
# t = 3.0536, df = 2.0005, p-value = 0.04628
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.08444013        Inf
# sample estimates:
#   mean of x mean of y 
# 1.9369473 0.0152417

t.test(df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN-mut3")], df1$betanin_ngmg[which(df1$sequence == "uninfiltrated")],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df1$betanin_ngmg[which(df1$sequence == "DODAa2-MN-mut3")] and df1$betanin_ngmg[which(df1$sequence == "uninfiltrated")]
# t = 4.1191, df = 4.0082, p-value = 0.007281
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.4446652       Inf
# sample estimates:
#   mean of x mean of y 
# 0.9363474 0.0152417 

# Yeast genomic integration, their sequence

df2 <- read.csv("Scerevisiae_GI.csv", header = T, stringsAsFactors = T)
str(df2)
df2$sequence <- factor(df2$sequence, c("DODAa1", "DODAa2-mut3", "DODAa2", "CYP76AD6"))

# plotting

means <- aggregate(fluorescence ~ sequence, data = df2, FUN = mean)
rownames(means) <- means[,1]
means <- means[,-1]
sds <- aggregate(fluorescence ~ sequence, data = df2, FUN = sd)
rownames(sds) <- sds[,1]
sds <- sds[,-1]

bp <- barplot(means, names.arg = c(expression(paste("DODA", alpha, "1")),
                                   expression(paste("DODA", alpha, "2-mut3")), 
                                   expression(paste("DODA", alpha, "2")),
                                   expression(paste("Control"))),
              col = "#FFCC00", border = "NA", yaxt = "n", ylim = c(0, 8200))

error.bar(bp, means, sds)
axis(2, c(0, 2000, 4000, 6000, 8000))

points(jitter(rep(.5, 4)), df2$fluorescence[ which(df2$sequence == "DODAa1") ], pch = 1)
points(jitter(rep(1.7, 4)), df2$fluorescence[ which(df2$sequence == "DODAa2-mut3") ], pch = 1)
points(jitter(rep(2.9, 4)), df2$fluorescence[ which(df2$sequence == "DODAa2") ], pch = 1)
points(jitter(rep(4.1, 4)), df2$fluorescence[ which(df2$sequence == "CYP76AD6") ], pch = 1)

title(main = expression(paste("DODA expression by genomic integration in ", italic("S. cerevisiae"))), adj = 0)
title(ylab = "Fluorescence")

X <- grconvertX(c(1.5, 5), from = "user", to = "ndc")
Y <- grconvertY(c(1000, 7000), from = "user", to = "ndc")

par(fig = c(X, Y), new = T)

bp <- barplot(means[-1], col = "#FFCC00", border = "NA", yaxt = "n", ylim = c(0, 400))
error.bar(bp, means[-1], sds[-1])
axis(2, at = c(0, 200, 400), labels = c(0, 200, 400))
points(jitter(rep(.5, 4)), df2$fluorescence[ which(df2$sequence == "DODAa2-mut3") ], pch = 1)
points(jitter(rep(1.7, 4)), df2$fluorescence[ which(df2$sequence == "DODAa2") ], pch = 1)
points(jitter(rep(3, 4)), df2$fluorescence[ which(df2$sequence == "CYP76AD6") ], pch = 1)


# stats

t.test(df2$fluorescence[ which(df2$sequence == "DODAa1") ], df2$fluorescence[ which(df2$sequence == "DODAa2") ],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df2$fluorescence[which(df2$sequence == "DODAa1")] and df2$fluorescence[which(df2$sequence == "DODAa2")]
# t = 94.602, df = 4.119, p-value = 2.48e-08
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   7464.781      Inf
# sample estimates:
#   mean of x mean of y 
# 7931.0442  295.6186 

t.test(df2$fluorescence[ which(df2$sequence == "DODAa1") ], df2$fluorescence[ which(df2$sequence == "DODAa2-mut3") ],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df2$fluorescence[which(df2$sequence == "DODAa1")] and df2$fluorescence[which(df2$sequence == "DODAa2-mut3")]
# t = 100.85, df = 3.2043, p-value = 5.062e-07
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   7404.5    Inf
# sample estimates:
#   mean of x mean of y 
# 7931.0442  354.2425

t.test(df2$fluorescence[ which(df2$sequence == "DODAa2-mut3") ], df2$fluorescence[ which(df2$sequence == "DODAa2") ],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df2$fluorescence[which(df2$sequence == "DODAa2-mut3")] and df2$fluorescence[which(df2$sequence == "DODAa2")]
# t = 1.6634, df = 4.0255, p-value = 0.08556
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   -16.37386       Inf
# sample estimates:
#   mean of x mean of y 
# 354.2425  295.6186 

t.test(df2$fluorescence[ which(df2$sequence == "DODAa2") ], df2$fluorescence[ which(df2$sequence == "CYP76AD6") ],
       var.equal = F, alternative = "g")

# Welch Two Sample t-test
# 
# data:  df2$fluorescence[which(df2$sequence == "DODAa2")] and df2$fluorescence[which(df2$sequence == "CYP76AD6")]
# t = 2.8352, df = 5.6062, p-value = 0.01603
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   35.43962      Inf
# sample estimates:
#   mean of x mean of y 
# 295.6186  179.7698

t.test(df2$fluorescence[ which(df2$sequence == "DODAa2-mut3") ], df2$fluorescence[ which(df2$sequence == "CYP76AD6") ],
       var.equal = F, alternative = "g")

# Welch Two Sample t-test
# 
# data:  df2$fluorescence[which(df2$sequence == "DODAa2-mut3")] and df2$fluorescence[which(df2$sequence == "CYP76AD6")]
# t = 6.1699, df = 4.6666, p-value = 0.001039
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   116.5792      Inf
# sample estimates:
#   mean of x mean of y 
# 354.2425  179.7698 

# Yeast high copy, their sequence

df3 <- read.csv("Scerevisiae_HC.csv", header = T, stringsAsFactors = T)
str(df3)
df3$sequence <- factor(df3$sequence, c("DODAa1", "DODAa2-mut3", "DODAa2", "pVV214"))

# plotting

means <- aggregate(fluorescence ~ sequence, data = df3, FUN = mean)
rownames(means) <- means[,1]
means <- means[,-1]
sds <- aggregate(fluorescence ~ sequence, data = df3, FUN = sd)
rownames(sds) <- sds[,1]
sds <- sds[,-1]

bp <- barplot(means, names.arg = c(expression(paste("DODA", alpha, "1")),
                                   expression(paste("DODA", alpha, "2-mut3")), 
                                   expression(paste("DODA", alpha, "2")),
                                   expression(paste("Control"))),
              col = "#FFCC00", border = "NA", yaxt = "n", ylim = c(0, 14500))

error.bar(bp, means, sds)
axis(2, at = c(0, 2000, 4000, 6000, 8000, 10000, 12000, 14000),
     labels = c(c(0, 2000, 4000, 6000, 8000, 10000, 12000, 14000)))
points(jitter(rep(.5, 4)), df3$fluorescence[ which(df3$sequence == "DODAa1") ], pch = 1)
points(jitter(rep(1.7, 4)), df3$fluorescence[ which(df3$sequence == "DODAa2-mut3") ], pch = 1)
points(jitter(rep(2.9, 4)), df3$fluorescence[ which(df3$sequence == "DODAa2") ], pch = 1)
points(jitter(rep(4.1, 4)), df3$fluorescence[ which(df3$sequence == "pVV214") ], pch = 1)

title(main = expression(paste("DODA expression by high-copy plasmid in ", italic("S. cerevisiae"))), adj = 0)
title(ylab = "Fluorescence")

X <- grconvertX(c(1.5, 5), from = "user", to = "ndc")
Y <- grconvertY(c(1500, 12000), from = "user", to = "ndc")

par(fig = c(X, Y), new = T)

bp <- barplot(means[-1], col = "#FFCC00", border = "NA", yaxt = "n", ylim = c(0, 1300))
error.bar(bp, means[-1], sds[-1])
axis(2, at = c(0, 400, 800, 1200))
points(jitter(rep(.5, 4)), df3$fluorescence[ which(df3$sequence == "DODAa2-mut3") ], pch = 1)
points(jitter(rep(1.7, 4)), df3$fluorescence[ which(df3$sequence == "DODAa2") ], pch = 1)
points(jitter(rep(3, 4)), df3$fluorescence[ which(df3$sequence == "pVV214") ], pch = 1)


# stats

t.test(df3$fluorescence[ which(df3$sequence == "DODAa1") ], df3$fluorescence[ which(df3$sequence == "DODAa2") ],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df3$fluorescence[which(df3$sequence == "DODAa1")] and df3$fluorescence[which(df3$sequence == "DODAa2")]
# t = 20.913, df = 3.0599, p-value = 0.0001052
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   10425.08      Inf
# sample estimates:
#   mean of x mean of y 
# 12758.948  1023.772 

t.test(df3$fluorescence[ which(df3$sequence == "DODAa1") ], df3$fluorescence[ which(df3$sequence == "DODAa2-mut3") ],
       var.equal = F, alternative = "greater")

# Welch Two Sample t-test
# 
# data:  df3$fluorescence[which(df3$sequence == "DODAa1")] and df3$fluorescence[which(df3$sequence == "DODAa2-mut3")]
# t = 20.933, df = 3.0448, p-value = 0.0001084
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   10420.39      Inf
# sample estimates:
#   mean of x mean of y 
# 12758.948  1027.511 

t.test(df3$fluorescence[ which(df3$sequence == "DODAa2-mut3") ], df3$fluorescence[ which(df3$sequence == "DODAa2") ],
       var.equal = F, alternative = "g")

# Welch Two Sample t-test
# 
# data:  df3$fluorescence[which(df3$sequence == "DODAa2-mut3")] and df3$fluorescence[which(df3$sequence == "DODAa2")]
# t = 0.050707, df = 5.8778, p-value = 0.4806
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   -140.0719       Inf
# sample estimates:
#   mean of x mean of y 
# 1027.511  1023.772 

t.test(df3$fluorescence[ which(df3$sequence == "DODAa2") ], df3$fluorescence[ which(df3$sequence == "pVV214") ],
       var.equal = F, alternative = "g")

# Welch Two Sample t-test
# 
# data:  df3$fluorescence[which(df3$sequence == "DODAa2")] and df3$fluorescence[which(df3$sequence == "pVV214")]
# t = 4.6861, df = 4.4519, p-value = 0.003617
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   321.9672      Inf
# sample estimates:
#   mean of x mean of y 
# 1023.7722  446.7893 

t.test(df3$fluorescence[ which(df3$sequence == "DODAa2-mut3") ], df3$fluorescence[ which(df3$sequence == "pVV214") ],
       var.equal = F, alternative = "g")

# Welch Two Sample t-test
# 
# data:  df3$fluorescence[which(df3$sequence == "DODAa2-mut3")] and df3$fluorescence[which(df3$sequence == "pVV214")]
# t = 4.8433, df = 4.1167, p-value = 0.003896
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   327.1785      Inf
# sample estimates:
#   mean of x mean of y 
# 1027.5109  446.7893 

