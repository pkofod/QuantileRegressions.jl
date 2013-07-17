# Author: Vincent Arel-Bundock (w/ changes by JMW)
# Contact: varel@umich.edu
# Date: 2013-07-13
# License: BSD-3
# Original code from the python statsmodels project
# https://github.com/statsmodels/statsmodels

# Inject dependency into caller
using DataFrames

module QuantileRegression
    # Need to export names you want caller to have access to.
    # Everything else is hidden to keep the Main namespace pure.
    export qreg

    # Access dependencies within module
    using DataFrames, Distributions

    # I like having types at the front, but that's just my taste
    type QRegModel
        beta::Vector{Float64}
        vcov::Matrix{Float64}
        mf::ModelFrame
    end

    # I like annotating function arguments to express meaning of parameters
    #  and to enforce correct function usage
    function bandwidth_hall_sheather(n::Integer, q::Real, alpha::Real)
        z = quantile(Normal(), q)
        num = 1.5 * pdf(Normal(), z)^2
        den = 2.0 * z^2 + 1
        h = n^(-1/3) * quantile(Normal(), (1. - alpha / 2.))^(2./3) * (num / den)^(1./3)
        return h
    end

    # Often better to start with a function that operates in-place
    #  then define a wrapper that allocates memory on each call
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

    # I really like using keyword arguments for defaults so that
    # function callers can modify them if need be.
    function qreg_coef(y::Vector, X::Matrix, q::Real = 0.5;
                       tol::Real = 1e-12, max_iter::Integer = 1_000,
                       threshold::Real = 1e-5)
        # What is this?
        n, p = size(X)

        # Allocate memory without spending time filling it
        beta0 = Array(Float64, p)
        beta = Array(Float64, p)

        # Keep a copy of xstar that we can reuse
        xstar = copy(X)

        # This ensures the loop runs at least once
        diff = Inf

        # Reuse this memory
        xtx = Array(Float64, p, p)
        xty = Array(Float64, p)
        xbeta = Array(Float64, n)
        resid = Array(Float64, n)

        for itr in 1:max_iter
            if diff > tol
                # Copy contents of beta into beta0 without allocating new space
                copy!(beta0, beta)

                # Compute without allocating memory
                # Note that this function name breaks Julia convention because
                # it should have an exclamation mark in it
                At_mul_B(xtx, xstar, X)
                # If you were using new memory, the idiomatic form would be
                # xtx = xstar' * X

                # Compute without allocating memory
                # Note that this function name breaks Julia convention because
                # it should have an exclamation mark in it
                At_mul_B(xty, xstar, y)
                # If you were using new memory, the idiomatic form would be
                # xty = xstar' * y

                # Instead of using pinv(), it's much more efficient to use
                # backslash. Sadly, most computer code that inverts a matrix
                # is numerically unstable.
                # The one thing that sucks is that I don't know how to do this
                # in place
                beta = xtx \ xty

                # In-place computation
                A_mul_B(xbeta, X, beta)
                for i in 1:n
                    @inbounds resid[i] = y[i] - xbeta[i]
                end
                # Slower, not in-place, idiomatic version
                # resid = y - X * beta

                # Changed to be in-place
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

                # In-place
                broadcast!(/, xstar, X, resid)
                # xstar = broadcast(/, X, resid)

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
