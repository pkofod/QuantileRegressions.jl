# Author: Vincent Arel-Bundock
# Contact: varel@umich.edu
# License: BSD-3
# Original code from the python statsmodels project
# https://github.com/statsmodels/statsmodels

function qreg_coef(y::Vector, X::Matrix, q::Real, s::IRLS)
    n, p = size(X)
    xstar = copy(X)
    diff = Inf

    beta0 = Array(Float64, p)
    beta = Array(Float64, p)
    xtx = Array(Float64, p, p)
    xty = Array(Float64, p)
    xbeta = Array(Float64, n)
    resid = Array(Float64, n)

    for itr in 1:s.maxIter
        if diff > s.tol
            copy!(beta0, beta)

            At_mul_B!(xtx, xstar, X)
            At_mul_B!(xty, xstar, y)
            beta = xtx \ xty
            A_mul_B!(xbeta, X, beta)

            for i in 1:n
                @inbounds resid[i] = y[i] - xbeta[i]
            end

            for i in 1:n
                if abs(resid[i]) < s.threshold
                    @inbounds resid[i] = sign(resid[i]) * s.threshold
                end
                if resid[i] < 0
                    @inbounds resid[i] = abs(q * resid[i])
                else
                    @inbounds resid[i] = abs((1 - q) * resid[i])
                end
            end

            broadcast!(/, xstar, X, resid)

            diff = norm(beta0-beta, Inf)
        end
    end
    return beta
end
