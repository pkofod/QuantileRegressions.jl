function npqreg(y, x, tau; m=50, h=2, xrange=nothing)
	if xrange = nothing
		xrange = collect(range(extrema(x)..., length=m))
	end
	z = copy(x)
	Z = hcat(fill(1, length(x)), z)
	
	ghat = copy(xrange)
	dghat = copy(xrange)
	for i in 1:length(xrange)
	  z .= x .- xrange[i]
	  Z[:, 2] .= z

	  w = pdf.(Normal(), z./h)

	  wy = w.*y
	  wZ = w.*Z 

	  r = qreg_coef(wy, wZ, tau, IRLS())

	  ghat[i] = r[1]
	  dghat[i] = r[2]
	end

	xrange, ghat, dghat
end
