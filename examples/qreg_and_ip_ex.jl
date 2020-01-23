using QuantileRegression
using DataFrames, CSV
cd(joinpath(dirname(pathof(QuantileRegression)), "..", "examples"))
df = CSV.read("engel.csv")


out_ip = qreg(foodexp~income, df, IP()) # or just qreg(foodexp~income, df)
out_irls = qreg(foodexp~income, df, IRLS())

println("...... Interior point ......")
println(out_ip)
println(".......... IRLS ............")
println(out_irls)
