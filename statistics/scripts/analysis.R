require(ggplot2)
require(multcomp)
require(sandwich)

#setwd("~/Dropbox/cary_projects/DODA/cary/bean_rebuttal/")
setwd("D:/Dropbox/cary_projects/DODA/cary/bean_rebuttal/bean_response/statistics/data/")
df1 <- read.csv("fig_4a.csv", header=T, stringsAsFactors = T)
str(df1)
df2 <- read.csv("fig_4b.csv", header=T, stringsAsFactors = T)
str(df2)
df3 <- read.csv("fig_4c.csv", header = T, stringsAsFactors = T)
str(df3)
df4 <- read.csv("fig_5.csv", header = T, stringsAsFactors = T)
str(df4)

ggplot(aes(x = sequence, y=betanin_equiv), data=df1) + 
  geom_bar(stat="identity")

t.test(df1$betanin_equiv[which(df1$sequence == "DODAa2-MN")], df1$betanin_equiv[which(df1$sequence == "DODAa2-MN-mut3")], var.equal = F)

# Welch Two Sample t-test
# 
# data:  df1$betanin_equiv[which(df1$sequence == "DODAa2-MN")] and df1$betanin_equiv[which(df1$sequence == "DODAa2-MN-mut3")]
# t = 1.5616, df = 2.2495, p-value = 0.2451
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.360461  3.196088
# sample estimates:
#   mean of x mean of y 
# 1.0912415 0.1734281 

ggplot(aes(x = sequence, y=log(fluorescence)), data=df2) + 
  geom_bar(stat="identity")

t.test(df2$fluorescence[which(df2$sequence == "DODAa2-MN")], df2$fluorescence[which(df2$sequence == "DODAa2-MN-mut3")], var.equal = F)

# Welch Two Sample t-test
# 
# data:  df2$fluorescence[which(df2$sequence == "DODAa2-MN")] and df2$fluorescence[which(df2$sequence == "DODAa2-MN-mut3")]
# t = 0.86118, df = 5.0017, p-value = 0.4285
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2035575  0.4086903
# sample estimates:
#   mean of x mean of y 
# 1.207716  1.105149 

ggplot(aes(x = sequence, y=log(fluorescence)), data=df3) + 
  geom_bar(stat="identity")

t.test(df3$fluorescence[which(df3$sequence == "DODAa2")], df3$fluorescence[which(df3$sequence == "DODAa2-mut3")], var.equal = F)

# Welch Two Sample t-test
# 
# data:  df3$fluorescence[which(df3$sequence == "DODAa2")] and df3$fluorescence[which(df3$sequence == "DODAa2-mut3")]
# t = -0.5698, df = 5.0237, p-value = 0.5933
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.259745  0.165378
# sample estimates:
#   mean of x mean of y 
# 1.179956  1.227139 

ggplot(aes(x = sequence, y=log(fluorescence)), data=df4) +
  geom_bar(stat="identity")

# name correspondence 
# yRG003  BvDODAa1
# yRG004  BvDODAa2
# yRG071	P1-BvDODAa1a2
# yRG069	P2-BvDODAa1a2
# yRG072	P3-BvDODAa1a2
# yRG070	P4-BvDODAa1a2
# yRG073	P5-BvDODAa1a2
# yRG053	P6-BvDODAa1a2
# yRG074	P7-BvDODAa1a2
# yRG054	P8-BvDODAa1a2
# yRG075	P9-BvDODAa1a2
# yRG055	P10-BvDODAa1a2
# yRG076	P11-BvDODAa1a2

# values are Fold change - construct/empty construct, and therefore must be log-transformed before comparison

fig5.fit1.log <- lm(log(fluorescence) ~ sequence, data=df) 
plot(fig5.fit1.log) # assumptions look okay
fig5.fit1.aov.log <- aov(log(fluorescence) ~ sequence, data=df)
summary(fig5.fit1.log)
summary(fig5.fit1.aov.log)
anova(fig5.fit1.log)

# Analysis of Variance Table
# 
# Response: log(fluorescence)
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# sequence  12 131.945 10.9954  957.98 < 2.2e-16 ***
#   Residuals 39   0.448  0.0115                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(fit1.aov.log)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = log(fluorescence) ~ sequence, data = df)
# 
# $sequence
# diff         lwr           upr     p adj
# yRG004-yRG003 -3.354839040 -3.62184824 -3.0878298362 0.0000000 *
# yRG053-yRG003  0.464716832  0.19770763  0.7317260361 0.0000235 *
# yRG054-yRG003  0.180177811 -0.08683139  0.4471870148 0.4787312
# yRG055-yRG003  0.571288534  0.30427933  0.8382977384 0.0000003 *
# yRG069-yRG003 -0.274684895 -0.54169410 -0.0076756912 0.0388291 *
# yRG070-yRG003 -0.266808632 -0.53381784  0.0002005717 0.0503270 
# yRG071-yRG003 -3.209968755 -3.47697796 -2.9429595507 0.0000000 *
# yRG072-yRG003 -2.892808248 -3.15981745 -2.6257990439 0.0000000 *
# yRG073-yRG003  0.436947715  0.16993851  0.7039569195 0.0000742 *
# yRG074-yRG003  0.581297021  0.31428782  0.8483062256 0.0000002 *
# yRG075-yRG003  0.776428853  0.50941965  1.0434380569 0.0000000 *
# yRG076-yRG003 -2.973006785 -3.24001599 -2.7059975813 0.0000000 *
# yRG053-yRG004  3.819555872  3.55254667  4.0865650765 0.0000000 *
# yRG054-yRG004  3.535016851  3.26800765  3.8020260552 0.0000000 *
# yRG055-yRG004  3.926127575  3.65911837  4.1931367788 0.0000000 *
# yRG069-yRG004  3.080154145  2.81314494  3.3471633492 0.0000000 *
# yRG070-yRG004  3.088030408  2.82102120  3.3550396121 0.0000000 *
# yRG071-yRG004  0.144870286 -0.12213892  0.4118794897 0.7784519 
# yRG072-yRG004  0.462030792  0.19502159  0.7290399965 0.0000262 *
# yRG073-yRG004  3.791786756  3.52477755  4.0587959599 0.0000000 *
# yRG074-yRG004  3.936136062  3.66912686  4.2031452660 0.0000000 *
# yRG075-yRG004  4.131267893  3.86425869  4.3982770973 0.0000000 *
# yRG076-yRG004  0.381832255  0.11482305  0.6488414591 0.0007039 *
# yRG054-yRG053 -0.284539021 -0.55154823 -0.0175298171 0.0278115 *
# yRG055-yRG053  0.106571702 -0.16043750  0.3735809065 0.9668167 
# yRG069-yRG053 -0.739401727 -1.00641093 -0.4723925231 0.0000000 *
# yRG070-yRG053 -0.731525464 -0.99853467 -0.4645162602 0.0000000 *
# yRG071-yRG053 -3.674685587 -3.94169479 -3.4076763826 0.0000000 *
# yRG072-yRG053 -3.357525080 -3.62453428 -3.0905158758 0.0000000 *
# yRG073-yRG053 -0.027769117 -0.29477832  0.2392400876 1.0000000 
# yRG074-yRG053  0.116580189 -0.15042901  0.3835893937 0.9375086 
# yRG075-yRG053  0.311712021  0.04470282  0.5787212250 0.0105836 *
# yRG076-yRG053 -3.437723617 -3.70473282 -3.1707144132 0.0000000 *
# yRG055-yRG054  0.391110724  0.12410152  0.6581199278 0.0004842 *
# yRG069-yRG054 -0.454862706 -0.72187191 -0.1878535017 0.0000353 *
# yRG070-yRG054 -0.446986443 -0.71399565 -0.1799772389 0.0000490 *
# yRG071-yRG054 -3.390146565 -3.65715577 -3.1231373613 0.0000000 *
# yRG072-yRG054 -3.072986059 -3.33999526 -2.8059768545 0.0000000 *
# yRG073-yRG054  0.256769905 -0.01023930  0.5237791089 0.0693221 
# yRG074-yRG054  0.401119211  0.13411001  0.6681284150 0.0003226 *
# yRG075-yRG054  0.596251042  0.32924184  0.8632602464 0.0000001 *
# yRG076-yRG054 -3.153184596 -3.42019380 -2.8861753919 0.0000000 *
# yRG069-yRG055 -0.845973430 -1.11298263 -0.5789642253 0.0000000 *
# yRG070-yRG055 -0.838097167 -1.10510637 -0.5710879625 0.0000000 *
# yRG071-yRG055 -3.781257289 -4.04826649 -3.5142480849 0.0000000 *
# yRG072-yRG055 -3.464096782 -3.73110599 -3.1970875781 0.0000000 *
# yRG073-yRG055 -0.134340819 -0.40135002  0.1326683853 0.8506608 
# yRG074-yRG055  0.010008487 -0.25700072  0.2770176914 1.0000000 
# yRG075-yRG055  0.205140319 -0.06186889  0.4721495227 0.2865888 
# yRG076-yRG055 -3.544295320 -3.81130452 -3.2772861155 0.0000000 *
# yRG070-yRG069  0.007876263 -0.25913294  0.2748854671 1.0000000 
# yRG071-yRG069 -2.935283860 -3.20229306 -2.6682746553 0.0000000 *
# yRG072-yRG069 -2.618123353 -2.88513256 -2.3511141486 0.0000000 *
# yRG073-yRG069  0.711632611  0.44462341  0.9786418149 0.0000000 *
# yRG074-yRG069  0.855981917  0.58897271  1.1229911210 0.0000000 *
# yRG075-yRG069  1.051113748  0.78410454  1.3181229523 0.0000000 *
# yRG076-yRG069 -2.698321890 -2.96533109 -2.4313126859 0.0000000 *
# yRG071-yRG070 -2.943160122 -3.21016933 -2.6761509182 0.0000000 *
# yRG072-yRG070 -2.625999616 -2.89300882 -2.3589904114 0.0000000 *
# yRG073-yRG070  0.703756348  0.43674714  0.9707655520 0.0000000 *
# yRG074-yRG070  0.848105654  0.58109645  1.1151148581 0.0000000 *
# yRG075-yRG070  1.043237485  0.77622828  1.3102466894 0.0000000 *
# yRG076-yRG070 -2.706198153 -2.97320736 -2.4391889488 0.0000000 *
# yRG072-yRG071  0.317160507  0.05015130  0.5841697110 0.0086587 *
# yRG073-yRG071  3.646916470  3.37990727  3.9139256744 0.0000000 *
# yRG074-yRG071  3.791265776  3.52425657  4.0582749805 0.0000000 *
# yRG075-yRG071  3.986397608  3.71938840  4.2534068118 0.0000000 *
# yRG076-yRG071  0.236961969 -0.03004723  0.5039711736 0.1254997 
# yRG073-yRG072  3.329755963  3.06274676  3.5967651676 0.0000000 *
# yRG074-yRG072  3.474105270  3.20709607  3.7411144737 0.0000000 *
# yRG075-yRG072  3.669237101  3.40222790  3.9362463051 0.0000000 *
# yRG076-yRG072 -0.080198537 -0.34720774  0.1868106668 0.9968559 
# yRG074-yRG073  0.144349306 -0.12265990  0.4113585103 0.7823380
# yRG075-yRG073  0.339481137  0.07247193  0.6064903416 0.0037287 *
# yRG076-yRG073 -3.409954501 -3.67696371 -3.1429452966 0.0000000 *
# yRG075-yRG074  0.195131831 -0.07187737  0.4621410356 0.3575910
# yRG076-yRG074 -3.554303807 -3.82131301 -3.2872946027 0.0000000 *
# yRG076-yRG075 -3.749435638 -4.01644484 -3.4824264340 0.0000000 *

# vs BvDODAa1

set.seed(12345)
summary(glht(fig5.fit1.aov.log, linfct=mcp(sequence="Dunnett")))

# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Dunnett Contrasts
# 
# 
# Fit: aov(formula = log(fluorescence) ~ sequence, data = df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
#   yRG004 - yRG003 == 0 -3.35484    0.07576 -44.285  < 0.001 ***
#   yRG053 - yRG003 == 0  0.46472    0.07576   6.134  < 0.001 ***
#   yRG054 - yRG003 == 0  0.18018    0.07576   2.378  0.15898    
#   yRG055 - yRG003 == 0  0.57129    0.07576   7.541  < 0.001 ***
#   yRG069 - yRG003 == 0 -0.27468    0.07576  -3.626  0.00779 ** 
#   yRG070 - yRG003 == 0 -0.26681    0.07576  -3.522  0.01057 *  
#   yRG071 - yRG003 == 0 -3.20997    0.07576 -42.373  < 0.001 ***
#   yRG072 - yRG003 == 0 -2.89281    0.07576 -38.186  < 0.001 ***
#   yRG073 - yRG003 == 0  0.43695    0.07576   5.768  < 0.001 ***
#   yRG074 - yRG003 == 0  0.58130    0.07576   7.673  < 0.001 ***
#   yRG075 - yRG003 == 0  0.77643    0.07576  10.249  < 0.001 ***
#   yRG076 - yRG003 == 0 -2.97301    0.07576 -39.245  < 0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

# vs BvDODAa2

df4$sequence <- relevel(df4$sequence, "yRG004")
set.seed(12345)
fig5.fit1.aov.log.bvda2.ref <- aov(log(fluorescence) ~ sequence, data=df4)
summary(glht(fig5.fit1.aov.log.bvda2.ref, linfct=mcp(sequence="Dunnett")))
df4$sequence <- relevel(df4$sequence, "yRG003") # reverse change

# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Dunnett Contrasts
# 
# 
# Fit: aov(formula = log(fluorescence) ~ sequence, data = df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
#   yRG003 - yRG004 == 0  3.35484    0.07576  44.285   <0.001 ***
#   yRG053 - yRG004 == 0  3.81956    0.07576  50.420   <0.001 ***
#   yRG054 - yRG004 == 0  3.53502    0.07576  46.664   <0.001 ***
#   yRG055 - yRG004 == 0  3.92613    0.07576  51.827   <0.001 ***
#   yRG069 - yRG004 == 0  3.08015    0.07576  40.659   <0.001 ***
#   yRG070 - yRG004 == 0  3.08803    0.07576  40.763   <0.001 ***
#   yRG071 - yRG004 == 0  0.14487    0.07576   1.912     0.37    
#   yRG072 - yRG004 == 0  0.46203    0.07576   6.099   <0.001 ***
#   yRG073 - yRG004 == 0  3.79179    0.07576  50.053   <0.001 ***
#   yRG074 - yRG004 == 0  3.93614    0.07576  51.959   <0.001 ***
#   yRG075 - yRG004 == 0  4.13127    0.07576  54.535   <0.001 ***
#   yRG076 - yRG004 == 0  0.38183    0.07576   5.040   <0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)
# 

# allow for heteroscedasticity

set.seed(12345)
summary(glht(fig5.fit1.aov.log, linfct=mcp(sequence="Dunnett"), vcov=vcovHC))

# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Dunnett Contrasts
# 
# 
# Fit: aov(formula = log(fluorescence) ~ sequence, data = df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
#   yRG004 - yRG003 == 0 -3.35484    0.08145 -41.189   <0.001 ***
#   yRG053 - yRG003 == 0  0.46472    0.09817   4.734   <0.001 ***
#   yRG054 - yRG003 == 0  0.18018    0.06399   2.816   0.0566 .  
#   yRG055 - yRG003 == 0  0.57129    0.07398   7.722   <0.001 ***
#   yRG069 - yRG003 == 0 -0.27468    0.10799  -2.544   0.1043    
#   yRG070 - yRG003 == 0 -0.26681    0.08705  -3.065   0.0314 *  
#   yRG071 - yRG003 == 0 -3.20997    0.06601 -48.630   <0.001 ***
#   yRG072 - yRG003 == 0 -2.89281    0.07892 -36.657   <0.001 ***
#   yRG073 - yRG003 == 0  0.43695    0.08935   4.890   <0.001 ***
#   yRG074 - yRG003 == 0  0.58130    0.09591   6.061   <0.001 ***
#   yRG075 - yRG003 == 0  0.77643    0.09749   7.964   <0.001 ***
#   yRG076 - yRG003 == 0 -2.97301    0.10199 -29.151   <0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

# vs BvDODAa2 

df4$sequence <- relevel(df4$sequence, "yRG004")
set.seed(12345)
fig5.fit1.aov.log.bvda2.ref <- aov(log(fluorescence) ~ sequence, data=df4)
summary(glht(fig5.fit1.aov.log.bvda2.ref, linfct=mcp(sequence="Dunnett"), vcov=vcovHC))
df4$sequence <- relevel(df4$sequence, "yRG003") # reverse change

# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Dunnett Contrasts
# 
# 
# Fit: aov(formula = log(fluorescence) ~ sequence, data = df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
#   yRG003 - yRG004 == 0  3.35484    0.08145  41.189  < 0.001 ***
#   yRG053 - yRG004 == 0  3.81956    0.09192  41.551  < 0.001 ***
#   yRG054 - yRG004 == 0  3.53502    0.05393  65.554  < 0.001 ***
#   yRG055 - yRG004 == 0  3.92613    0.06547  59.968  < 0.001 ***
#   yRG069 - yRG004 == 0  3.08015    0.10235  30.096  < 0.001 ***
#   yRG070 - yRG004 == 0  3.08803    0.07994  38.629  < 0.001 ***
#   yRG071 - yRG004 == 0  0.14487    0.05630   2.573  0.10907    
#   yRG072 - yRG004 == 0  0.46203    0.07100   6.508  < 0.001 ***
#   yRG073 - yRG004 == 0  3.79179    0.08244  45.997  < 0.001 ***
#   yRG074 - yRG004 == 0  3.93614    0.08951  43.974  < 0.001 ***
#   yRG075 - yRG004 == 0  4.13127    0.09120  45.299  < 0.001 ***
#   yRG076 - yRG004 == 0  0.38183    0.09599   3.978  0.00297 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)
