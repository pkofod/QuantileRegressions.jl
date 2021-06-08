using GLM, QuantileRegressions, CSV, Winston, StatsBase, DataFrames

# Load data
url = "http://vincentarelbundock.github.io/Rdatasets/csv/quantreg/engel.csv"
df = CSV.read(download(url), DataFrame)
ResultQR = qreg(@formula(foodexp ~ income), df, .5)

# Fit quantile regression for a bunch of different quantiles
QNum = 35; # Number of equally spaced quantiles
DataPlot = reduce(vcat, [[coeftable(qreg(@formula(foodexp ~ income), df, i/(QNum+1), IP())).cols[2][2] coeftable(qreg(@formula(foodexp ~ income), df, i/(QNum+1), IP())).cols[3][2]] for i in 1:QNum])

# Fit OLS model to compare
ResultLM = lm(@formula(foodexp ~ income), df) # Fit the model using OLS from the GLM-package
ols = vcat([[coeftable(ResultLM).cols[1][2] coeftable(ResultLM).cols[2][2]] for i in 1:QNum]...)

PlotDF = DataFrame(
    X = range(1, length=QNum, stop=QNum)/(QNum+1),
    Y = DataPlot[:,1],
    ols = ols[:,1],
    Ymin = DataPlot[:,1] - 1.96 * DataPlot[:,2],
    Ymax = DataPlot[:,1] + 1.96 * DataPlot[:,2],
    YminLM = ols[:,1] - 1.96 * ols[:,2],
    YmaxLM = ols[:,1] + 1.96 * ols[:,2]
)

x = PlotDF[!, :X]
p = FramedPlot()
βτ   = Curve(x, PlotDF[!, :Y])
βols = Curve(x, PlotDF[!, :ols] , "color", "red")
ci_low  = Curve(x, PlotDF[!, :Ymin], "type" , "dash")
ci_high = Curve(x, PlotDF[!, :Ymax], "type" , "dash")
ci_lowLM  = Curve(x, PlotDF[!, :YminLM], "type" , "dash", "color", "red")
ci_highLM = Curve(x, PlotDF[!, :YmaxLM], "type" , "dash", "color", "red")
setattr(βτ,"label","β(τ)")
setattr(βols, "label", "β(ols)")
lgnd = Legend(.8,.2,[βτ, βols])
add(p, βτ, βols, ci_low, ci_high, ci_lowLM, ci_highLM, lgnd)
setattr(p, "title" , "Quantile regression: Food Expenditure ~ Income")
setattr(p, "xlabel", "Quantiles")
setattr(p, "ylabel", "Coefficient on Income")

# savefig(p, Pkg.dir("QuantileRegressions")"/examples/qreg_example_plot.png") # not working
