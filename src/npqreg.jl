function npqreg(y, x, tau, method=IP(); m=50, h=2, xrange=nothing)
    T = eltype(x)
    if xrange === nothing
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
        widx = findall(x-> x > eps(T), w)
        if length(widx) > 1
            w = w[widx]
            wy = w.*y[widx]
            wZ = w.*Z[widx, :]
            r = qreg_coef(wy, wZ, tau, method)
            ghat[i] = r[1]
            dghat[i] = r[2]
        else
            ghat[i] = y[widx[1]]
            dghat[i] = NaN
        end
    end

    xrange, ghat, dghat
end
