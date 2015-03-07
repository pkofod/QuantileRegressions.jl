using GLM, QuantileRegression, Gadfly

# Load data
url = "http://vincentarelbundock.github.io/Rdatasets/csv/quantreg/engel.csv"
# TODO: Make automatic URL downloads work
dat = readtable("engel.csv")

# Fit least absolute deviation model (quantile  = .5)
res = qreg(foodexp~income, dat, .5)
coefs = coeftable(res)
show(coefs)
# Fit quantile regression for a bunch of different quantiles
data_plot = reduce(vcat,[coeftable(qreg(foodexp ~ income, dat, i/20)).mat[2,1:2] for i in 1:19])

# Fit OLS model to compare
res_lm = lm(foodexp~income, dat)
ols = reduce(vcat,rep(coeftable(res_lm).mat[2,1:2],19))

Ydata = [ols[:,1];data_plot[:,1]]
Yse   = [ols[:,2];data_plot[:,2]]

PlotDataFrame =  DataFrame(X = [linspace(1,19,19)/20;linspace(1,19,19)/20],
                         Y = Ydata,
                         Ymin = Ydata - 1.96 * Yse,
                         Ymax = Ydata + 1.96 * Yse, 
                         Estimator = [["OLS" for i = 1:19];["Q(0.5)" for i = 1:19]])

plot(layer(PlotDataFrame, x = "X", y = "Ymin",Geom.line,color = "Estimator"),
     layer(PlotDataFrame, x = "X", y = "Ymax",Geom.line,color = "Estimator"),
     layer(PlotDataFrame, x = "X", y = "Y", color = "Estimator",Geom.line), 
     Guide.XLabel("Income"), Guide.YLabel("Food Expenditure"))
