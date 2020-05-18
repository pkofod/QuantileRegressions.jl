using QuantileRegressions, CSV
@testset "npqreg" begin
	mycyle = CSV.read("../examples/mcycle.csv")
	x = mycyle.times
	y = mycyle.accel

	for t = range(0.01,0.99,length=99)
		res = QuantileRegressions.npqreg(Array(y), Array(x), t; h=2)
		xrange = sort(unique(x))
		res = QuantileRegressions.npqreg(Array(y), Array(x), t; xrange=xrange, h=0.5)
	end
	res = QuantileRegressions.npqreg([y[1]], [x[1]], 0.4; xrange=[x[1]], h=2)
	res = QuantileRegressions.npqreg([y[1]], [x[1]], 0.8; xrange=[x[1]], h=0.1)
end
