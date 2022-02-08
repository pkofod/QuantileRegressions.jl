Data = CSV.read("engel.csv", DataFrame)

out_ip_true = qreg(@formula(foodexp~income), Data, IP(true))
out_ip = qreg(@formula(foodexp~income), Data, IP())
out_ip_false = qreg(@formula(foodexp~income), Data, IP(false))
out_irls = qreg(@formula(foodexp~income), Data, IRLS())

@test norm( out_ip_true.model.beta - [81.4822; 0.560181]) < 1e-4
@test norm( out_ip.model.beta - [81.4822; 0.560181]) < 1e-4
@test norm( out_ip_false.model.beta - [81.4822; 0.560181]) < 1e-4
@test norm( out_irls.model.beta - [81.4823; 0.560181]) < 1e-4

#=
...... Interior point ......
TableRegressionModel{QRegModel,Float64}:

Coefficients:
             Estimate Std.Error t value
(Intercept)   81.4822   14.6345 5.56783
income       0.560181 0.0131756 42.5164

.......... IRLS ............
TableRegressionModel{QRegModel,Float64}:

Coefficients:
             Estimate Std.Error t value
(Intercept)   81.4823   14.6345 5.56783
income       0.560181 0.0131756 42.5164
=#
