using QuantileRegressions, Distributions, CSV
mycyle = CSV.read("mcycle.csv")
x = mycyle.times
y = mycyle.accel

for t = range(0.01,0.99,length=99)
res = QuantileRegressions.npqreg(Array(y), Array(x), t; h=2)
plot!(res[1], res[2]; label="tau=$t")
end