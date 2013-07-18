# Author: Vincent Arel-Bundock (w/ changes by JMW)
# Contact: varel@umich.edu
# License: BSD-3
# Original code from the python statsmodels project
# https://github.com/statsmodels/statsmodels

using DataFrames

module QuantileRegression
    export qreg

    using DataFrames, Distributions

    type QRegModel
        beta::Vector{Float64}
        vcov::Matrix{Float64}
        mf::ModelFrame
    end

    function bandwidth_hall_sheather(n::Integer, q::Real, alpha::Real)
        z = quantile(Normal(), q)
        num = 1.5 * pdf(Normal(), z)^2
        den = 2.0 * z^2 + 1
        h = n^(-1/3) * quantile(Normal(), (1. - alpha / 2.))^(2./3) * (num / den)^(1./3)
        return h
    end

    # operate in-place to save on memory allocation
    function kernel_epanechnikov!(u::Vector, x::Vector)
        for i in 1:length(x)
            if abs(x[i]) <= 1
                # The @inbounds macro tells the compiler it can remove
                # bounds-checking for array indexing
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

    function maxdiff(x::Array, y::Array)
        m = typemin(Float64)
        for i in length(x)
            t = abs(x[i] - y[i])
            if t > m
                m = t
            end
        end
        return m
    end

    function qreg_coef(y::Vector, X::Matrix, q::Real = 0.5;
                       tol::Real = 1e-12, max_iter::Integer = 1_000,
                       threshold::Real = 1e-5)
        n, p = size(X)
        xstar = copy(X)
        diff = Inf

        beta0 = Array(Float64, p)
        beta = Array(Float64, p)
        xtx = Array(Float64, p, p)
        xty = Array(Float64, p)
        xbeta = Array(Float64, n)
        resid = Array(Float64, n)

        for itr in 1:max_iter
            if diff > tol
                copy!(beta0, beta)

                At_mul_B(xtx, xstar, X)
                At_mul_B(xty, xstar, y)
                beta = xtx \ xty
                A_mul_B(xbeta, X, beta)

                for i in 1:n
                    @inbounds resid[i] = y[i] - xbeta[i]
                end

                for i in 1:n
                    if abs(resid[i]) < threshold
                        @inbounds resid[i] = sign(resid[i]) * threshold
                    end
                    if resid[i] < 0
                        @inbounds resid[i] = abs(q * resid[i])
                    else
                        @inbounds resid[i] = abs((1 - q) * resid[i])
                    end
                end

                broadcast!(/, xstar, X, resid)

                diff = maxdiff(beta0, beta)
            end
        end
        return beta
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
        # This pinv() makes me uncomfortable
        xtxi = pinv(X' * X)
        xtdx = broadcast(*, X, d)' * X
        vcov = (xtxi * xtdx) * xtxi
        return vcov
    end

    function qreg(f::Expr, df::AbstractDataFrame, q::Real = 0.5)
        mf = ModelFrame(f, df)
        mm = ModelMatrix(mf)
        mr = model_response(mf)
        coef = qreg_coef(mr, mm.m, q)
        vcov = qreg_vcov(mr, mm.m, coef, q)
        return QRegModel(coef, vcov, mf)
    end

    coef(x::QRegModel) = x.beta

    vcov(x::QRegModel) = x.vcov

    function DataFrames.describe(x::QRegModel)
        # TODO: Insert column names
        se = sqrt(diag(x.vcov))
        out = DataFrame()
        out["Estimate"] = x.beta
        out["Std.Error"] = se
        out["t value"] = x.beta ./ se
        out["2.5%"] = x.beta - 1.96 * se
        out["97.5%"] = x.beta + 1.96 * se
        return out
    end
end
