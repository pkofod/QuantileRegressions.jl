using GLM, Winston, QReg

# Load data
url = "http://vincentarelbundock.github.io/Rdatasets/csv/quantreg/engel.csv"
dat = readtable("/Users/USERNAME/Downloads/engel.csv")

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
