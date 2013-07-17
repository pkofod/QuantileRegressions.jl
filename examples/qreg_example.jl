using GLM, Winston, QReg

# Load data
url = "http://vincentarelbundock.github.io/Rdatasets/csv/quantreg/engel.csv"
dat = readtable("engel.csv")

# Fit least absolute deviation model (quantile =.5)
res = qreg(:(foodexp~income), dat, .5)
summary(res)

# Fit quantile regression for a bunch of different quantiles
dat_plot = [summary(qreg(:(foodexp ~ income), dat, i/10))[2,:] for i in 1:9]
dat_plot = reduce(vcat, dat_plot)

# Fit OLS model to compare
res_lm = lm(:(foodexp~income), dat)
ols = coeftable(res_lm)[2,1]
ols = rep(ols, 9)

# Plot results
x = [i/10 for i=1:9]
y = dat_plot["Estimate"]
p = FramedPlot()
add(p, Curve(x, dat_plot["Estimate"]))
add(p, Curve(x, dat_plot["2.5%"], "type", "dash"))
add(p, Curve(x, dat_plot["97.5%"], "type", "dash"))
add(p, Curve(x, ols, "color", "red"))
setattr(p, "title", "Quantile regression: Food Expenditure ~ Income\n Red bar: OLS coefficient")
setattr(p, "xlabel", "Quantiles")
setattr(p, "ylabel", "Coefficient on Income")
file(p, "/Users/USERNAME/qreg_example_plot.png")

# > summary(res)
# 2x5 DataFrame:
#         Estimate Std.Error t value     2.5%    97.5%
# [1,]     81.4823   14.6345 5.56783  52.7987  110.166
# [2,]    0.560181 0.0131756 42.5164 0.534356 0.586005
# insheet using engel.csv
# qreg foodexp income, vce(iid, kernel(epan2))
 
# . qreg foodexp income, vce(iid, kernel(epan2))
# Iteration  1:  WLS sum of weighted deviations =  18051.196
# 
# Iteration  1: sum of abs. weighted deviations =  18035.127
# Iteration  2: sum of abs. weighted deviations =  17815.249
# Iteration  3: sum of abs. weighted deviations =  17567.085
# Iteration  4: sum of abs. weighted deviations =  17559.932
# 
# Median regression                                    Number of obs =       235
#   Raw sum of deviations 46278.06 (about 582.54126)
#   Min sum of deviations 17559.93                     Pseudo R2     =    0.6206
# 
# ------------------------------------------------------------------------------
#      foodexp |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
# -------------+----------------------------------------------------------------
#       income |   .5601805   .0131763    42.51   0.000     .5342206    .5861403
#        _cons |   81.48233   14.63518     5.57   0.000     52.64815    110.3165
# ------------------------------------------------------------------------------
