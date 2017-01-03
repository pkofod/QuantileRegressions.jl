Data = readtable("engel.csv")

out_ip = qreg(foodexp~income, Data, IP())
out_irls = qreg(foodexp~income, Data, IRLS())

@assert norm( out_ip.model.beta - [81.4822; 0.560181]) < 1e-4
@assert norm( out_irls.model.beta - [81.4823; 0.560181]) < 1e-4

#=
...... Interior point ......
DataFrameRegressionModel{QRegModel,Float64}:

Coefficients:
             Estimate Std.Error t value
(Intercept)   81.4822   14.6345 5.56783
income       0.560181 0.0131756 42.5164

.......... IRLS ............
DataFrameRegressionModel{QRegModel,Float64}:

Coefficients:
             Estimate Std.Error t value
(Intercept)   81.4823   14.6345 5.56783
income       0.560181 0.0131756 42.5164
=#
