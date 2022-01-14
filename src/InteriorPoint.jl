# A Frisch-Newton algorithm as described in http://projecteuclid.org/euclid.ss/1030037960
# qreg_ip_coef is a modified lp_fnm of Daniel Morillo & Roger Koenker
# Translated from Ox to Matlab by Paul Eilers 1999
# Modified by Roger Koenker 2000
# More changes by Paul Eilers 2004
# QuantileRegressions.jl developers 2015
#
# The functions in this file is released under the BSD-3 license, although the original
# versions were not. Permissing was given to Patrick Kofod Mogensen by Roger Koenker and
# Paul Eilers on April 24 2015 via E-mail. Original mail correspondences can be provided
# upon request.

struct IP <: Solver
    cholesky::Bool
end
IP() = IP(false)

function bound!(d, x, dx)
# Fill vector with allowed step lengths
# Replace with -x/dx for negative dx
    d .= 1e20
    @inbounds for i = 1:length(dx)
        if dx[i] < 0.0
            d[i] = -x[i] / dx[i]
        end
    end
    return d
end

function qreg_coef(Y, X::Matrix, p, method::IP)
# input: X is an n x k matrix of exogenous regressors,
#        Y is an n x 1 vector of outcome variables
#        p \in (0,1) is the quantile of interest
# Output: p^{th} regression quantiles.
# Construct the dual problem of quantile regression
c = -Y
T = eltype(X)

m = size(X, 1)

x = fill((1 - p), m)
b = X'x
# Solve a linear program by the interior point method:
# min(c * u), s.t. A * x = b and 0 < x < u
# An initial feasible solution has to be provided as x

# Set some constants
beta = 0.99995
small = sqrt(eps(eltype(c)))
max_it = 500
n, m = size(X)
# Generate inital feasible point
s = 1 .- x
if method.cholesky
    y = -((X'X)\(X'Y))
else
    y = -X\Y
end
dy = copy(y)
r = c - X*y
# z and w are the "residuals" but on either
# side of the "check" we need to flip the sign
#    BLAS.axpy!(small, (r .== 0.0).*1.0, r)
#    z = r .* (r .> small)
#    w = abs.(z .- r)
z = copy(r)
w = copy(r)
for i ∈ eachindex(r)
    ri = r[i]
    z[i] = max(ri, 0)
    w[i] = max(-ri, 0)
    if abs(ri) <= small
        z[i] = z[i] + small
        w[i] = w[i] + small
    end
end

gap = dot(z, x) + dot(s, w)

# set up caches
Xtmp = copy(X)
if method.cholesky
    XtX  = similar(X, m, m)
end
xinv, xi = copy(x), copy(x)
sinv = copy(s)
q = copy(z)
dx, ds, dz, dw = copy(w), copy(w), copy(w), copy(w)
fx, fs, fz, fw = copy(w), copy(w), copy(w), copy(w)
dxdz, dsdw = copy(w), copy(w)
tmp = copy(w)
if method.cholesky
    rtmp = copy(r)
    Xtqr = similar(r, m)
else
    qr_tmp = copy(q)
    rhs = copy(r)
end
# Start iterations
for it = 1:max_it
    #   Compute affine step
    @. q = 1 / (z / x + w / s)
    @. r = z - w

    if method.cholesky
        @. Xtmp = q * X
        mul!(XtX, Xtmp', X)
        F = cholesky!(Symmetric(XtX))
        mul!(Xtqr, Xtmp', r)
        dy .= F\Xtqr
        #ldiv!(dy, F, Xtqr)
    else
        @. qr_tmp = sqrt(q)
        Q = Diagonal(qr_tmp) # Very efficient to do since Q diagonal
        AQtF = qr(mul!(Xtmp, Q, X)) # PE 2004
        @. rhs = qr_tmp*r
        dy .= AQtF\rhs
    end

    mul!(tmp, X, dy)
    @. dx = q*(tmp - r)
    @. ds = -dx
    @. dz = -z * (1 + dx / x)
    @. dw = -w * (1 + ds / s)
    # Compute maximum allowable step lengths
    bound!(fx, x, dx)
    bound!(fs, s, ds)
    bound!(fw, w, dw)
    bound!(fz, z, dz)
    @. tmp = min(fx, fs)
    fp = min(beta*minimum(tmp), 1)
    @. tmp = min(fw, fz)
    fd = min(beta*minimum(tmp), 1)
    # If full step is feasible, take it. Otherwise modify it
    if min(fp, fd) < 1.0
        # Update mu
        g0 = dot(z, x) + dot(w, s)
        gγPγD = T(0)
        for i = 1:length(z)
          @inbounds gγPγD += (z[i] + fd*dz[i])*(x[i] + fp*dx[i]) + (w[i] + fd*dw[i])*(s[i]+fp*ds[i])
        end
        mu = (gγPγD / g0)^3 * g0/(2n)
        # Compute modified step
        @. dxdz = dx * dz
        @. dsdw = ds * dw
        @. xinv = 1 / x
        @. sinv = 1 / s
        @. xi = mu * (xinv - sinv)
        if method.cholesky
            @. rtmp = r + dxdz - dsdw - xi
            mul!(Xtqr, Xtmp', rtmp)
            dy .= F\Xtqr
        else
            @. qr_tmp = qr_tmp * (dxdz .- dsdw .- xi)
            BLAS.axpy!(1.0, qr_tmp, rhs) # no gemv-wrapper gemv(Q, (dxdz - dsdw - xi), rhs,1,1,n)?
            dy .= AQtF\rhs
        end

        mul!(tmp, X, dy)
        @. dx = q * (tmp + xi - r - dxdz + dsdw)
        @. ds = -dx
        @. dz = (mu - z * dx)*xinv - z - dxdz
        @. dw = (mu - w * ds)*sinv - w - dsdw

        # Compute maximum allowable step lengths
        bound!(fx, x, dx)
        bound!(fs, s, ds)
        bound!(fw, w, dw)
        bound!(fz, z, dz)
        @. tmp = min(fx, fs)
        fp = min(beta*minimum(tmp), 1)
        @. tmp = min(fw, fz)
        fd = min(beta*minimum(tmp), 1)
    end

    # Take the steps
    BLAS.axpy!(fp, dx, x)
    BLAS.axpy!(fp, ds, s)
    BLAS.axpy!(fd, dy, y)
    BLAS.axpy!(fd, dw, w)
    BLAS.axpy!(fd, dz, z)

    gap = dot(z, x) + dot(s, w)

    if gap < small || isnan(gap)
        break
    end
end

return -y
end
