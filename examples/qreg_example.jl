using GLM, QuantileRegression, Winston

# Load data
url = "http://vincentarelbundock.github.io/Rdatasets/csv/quantreg/engel.csv"
Data = readtable("engel.csv")

# Fit least absolute deviation model (quantile  = .5)
ResultQR = qreg(foodexp~income, Data, .5)
β = coeftable(ResultQR)
show(β)

# Fit quantile regression for a bunch of different quantiles
QNum = 35; # Number of equally spaced quantiles
DataPlot = reduce(vcat,[coeftable(qreg(foodexp ~ income, Data, i/(QNum+1))).mat[2,1:2] for i in 1:QNum])

# Fit OLS model to compare
ResultLM = lm(foodexp~income, Data) # Fit the model using OLS from the GLM-package
ols = reduce(vcat,rep(coeftable(ResultLM).mat[2,1:2],QNum)) 

PlotDF =  DataFrame(X = linspace(1,QNum,QNum)/(QNum+1),
                         Y = DataPlot[:,1],
                         ols = ols[:,1],
                         Ymin = DataPlot[:,1] - 1.96 * DataPlot[:,2],
                         Ymax = DataPlot[:,1] + 1.96 * DataPlot[:,2])

x = PlotDF[:X]
p = FramedPlot()
βτ   = Curve(x, PlotDF[:Y])
βols = Curve(x, PlotDF[:ols] , "color", "red")
ci_low  = Curve(x, PlotDF[:Ymin], "type" , "dash")
ci_high = Curve(x, PlotDF[:Ymax], "type" , "dash")
setattr(βτ,"label","β(τ)")
setattr(βols, "label", "β(ols)")
lgnd = Legend(.8,.2,{βτ, βols})
add(p, βτ, βols, ci_low, ci_high, lgnd)
setattr(p, "title" , "Quantile regression: Food Expenditure ~ Income")
setattr(p, "xlabel", "Quantiles")
setattr(p, "ylabel", "Coefficient on Income")
savefig(p, "./qreg_example_plot.png")
