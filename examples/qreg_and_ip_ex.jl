using QuantileRegression
using DataFrames
Data = readtable("engel.csv")


out_ip = qreg(foodexp~income, Data, method = :ip) # or just qreg(foodexp~income, Data)
out_irls = qreg(foodexp~income, Data, method = :irls)

println("...... Interior point ......")
println(out_ip)
println(".......... IRLS ............")
println(out_irls)
