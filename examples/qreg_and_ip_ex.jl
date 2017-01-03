using QuantileRegression
using DataFrames
cd(Pkg.dir("QuantileRegression")*"/examples")
df = readtable("engel.csv")


out_ip = qreg(foodexp~income, df, IP()) # or just qreg(foodexp~income, df)
out_irls = qreg(foodexp~income, df, IRLS())

println("...... Interior point ......")
println(out_ip)
println(".......... IRLS ............")
println(out_irls)
