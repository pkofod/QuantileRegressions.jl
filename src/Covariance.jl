# Author: Vincent Arel-Bundock
# Contact: varel@umich.edu
# License: BSD-3
# Original code from the python statsmodels project
# https://github.com/statsmodels/statsmodels

function bandwidth_hall_sheather(n::Integer, q::Real, alpha::Real)
    z = quantile(Normal(), q)
    num = 1.5 * pdf(Normal(), z)^2
    den = 2.0 * z^2 + 1
    h = n^(-1/3) * quantile(Normal(), (1 - alpha / 2))^(2/3) * (num / den)^(1/3)
    return h
end

function kernel_epanechnikov!(u::Vector, x::Vector)
    for i in 1:length(x)
        if abs(x[i]) <= 1
            @inbounds u[i] = 0.75 * (1 - x[i]^2)
        else
            @inbounds u[i] = 0
        end
    end
    return
end

function kernel_epanechnikov(x::Vector)
u = similar(x)
kernel_epanechnikov!(u, x)
return u
end

function qreg_vcov(y::Vector, X::Matrix, beta::Vector, q::Real)
    resid = y - X * beta
    n = length(resid)
    iqre = quantile(resid, .75) - quantile(resid, .25)
    h = bandwidth_hall_sheather(n, q, .05)
    h1 = min(std(y), iqre / 1.34)
    h2 = quantile(Normal(), q + h) - quantile(Normal(), q - h)
    h = h1 * h2
    # TODO: This line could be optimized, but is a hassle
    fhat0 = 1 / (n * h) * sum(kernel_epanechnikov(resid / h))
    d = copy(resid)
    for i in 1:n
        if d[i] > 0
            @inbounds d[i] = (q / fhat0)^2
        else
            @inbounds d[i] = ((1 - q) / fhat0)^2
        end
    end

    xtx = X' * X
    xtdx = broadcast(*, X, d)' * X
    vcov = (xtx \ xtdx) / xtx

    return vcov
end
