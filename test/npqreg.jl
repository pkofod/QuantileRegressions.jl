using QuantileRegressions, CSV
@testset "npqreg" begin
	mycyle = CSV.read("../examples/mcycle.csv", DataFrame)
	x = mycyle.times
	y = mycyle.accel

	for t = range(0.01,0.99,length=99)
		res = QuantileRegressions.npqreg(Array(y), Array(x), t; h=2)
		@test all(!isnan, res[2])
		@test all(!isnan, res[3])
		xrange = sort(unique(x))
		res = QuantileRegressions.npqreg(Array(y), Array(x), t; xrange=xrange, h=0.5)
	end
	res = QuantileRegressions.npqreg([y[1]], [x[1]], 0.4; xrange=[x[1]], h=2)
	@test !isnan(res[2][1])
	@test isnan(res[3][1])
	res = QuantileRegressions.npqreg([y[1]], [x[1]], 0.8; xrange=[x[1]], h=0.0001)
	@test !isnan(res[2][1])
	@test isnan(res[3][1])
	res = QuantileRegressions.npqreg([y[1]], [x[1]], 0.8; xrange=[x[1].+1], h=0.1)
	@test isnan(res[2][1])
	@test isnan(res[3][1])

	@test all(QuantileRegressions.npqreg([1.0,2.0,4.0,5.0], [0.0,0.0,6.0,12.0], 0.3; h=0.001, xrange=[0.0,6.0,12.0])[2] .≈ [1.0, 4.0, 5.0])
end
