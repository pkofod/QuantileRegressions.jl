# Author: Vincent Arel-Bundock
# Contact: varel@umich.edu
# License: BSD-3
# Original code from the python statsmodels project
# https://github.com/statsmodels/statsmodels

function qreg_coef(y, X::Matrix, q::Real, method::IRLS)
    n, p = size(X)
    xstar = copy(X)
    diff = Inf

    beta0 = zeros(p)
    beta  = zeros(p)
    xtx   = Array{Float64}(undef, p, p)
    xty   = Array{Float64}(undef, p)
    xbeta = Array{Float64}(undef, n)
    resid = Array{Float64}(undef, n)

    for itr in 1:method.maxiter
        if diff > method.tol
            copyto!(beta0, beta)

            mul!(xtx, xstar', X)
            mul!(xty, xstar', y)
            beta .= xtx \ xty
            mul!(xbeta, X, beta)

            @. resid = y - xbeta

            for i in 1:n
                if abs(resid[i]) < method.threshold
                    @inbounds resid[i] = sign(resid[i]) * method.threshold
                end
                if resid[i] < 0
                    @inbounds resid[i] = abs(q * resid[i])
                else
                    @inbounds resid[i] = abs((1 - q) * resid[i])
                end
            end

            xstar .= X ./ resid

            diff = norm(beta0 - beta, Inf)
        end
    end
    return beta
end
