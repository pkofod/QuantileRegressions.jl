using DataFrames

module QuantileRegression

    import StatsBase.coeftable, DataFrames.DataFrameRegressionModel
    export qreg, coef, vcov, stderr

    using DataFrames, Distributions

    type QRegModel
        beta::Vector{Float64}
        vcov::Matrix{Float64}
        stderr::Vector{Float64}
    end

    function bandwidth_hall_sheather(n::Integer, q::Real, alpha::Real)
        z = quantile(Normal(), q)
        num = 1.5 * pdf(Normal(), z)^2
        den = 2.0 * z^2 + 1
        h = n^(-1/3) * quantile(Normal(), (1. - alpha / 2.))^(2./3) * (num / den)^(1./3)
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

    function qreg_coef(y::Vector, X::Matrix, q::Real = 0.5;
                       tol::Real = 1e-12, maxIter::Integer = 1_000,
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

        for itr in 1:maxIter
            if diff > tol
                copy!(beta0, beta)

                At_mul_B!(xtx, xstar, X)
                At_mul_B!(xty, xstar, y)
                beta = xtx \ xty
                A_mul_B!(xbeta, X, beta)

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

                diff = norm(beta0-beta, Inf)
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
    
        xtx = X' * X
        xtdx = broadcast(*, X, d)' * X
        vcov = (xtx \ xtdx) / xtx
        
        return vcov
    end

    function qreg(f::Formula, df::AbstractDataFrame, q::Real = 0.5)
        mf = ModelFrame(f, df)
        mm = ModelMatrix(mf)
        mr = model_response(mf)
        coef = qreg_coef(mr, mm.m, q)
        vcov = qreg_vcov(mr, mm.m, coef, q)
        stderr = sqrt(diag(vcov))
        return DataFrameRegressionModel(QRegModel(coef, vcov, stderr),mf,mm)
    end

    coef(x::QRegModel) = x.beta

    vcov(x::QRegModel) = x.vcov

    stderr(x::QRegModel) = x.stderr

    function coeftable (mm::QRegModel)
        cc = coef(mm)
        se = stderr(mm)
        tt = cc./se
        CoefTable(hcat(cc,se,tt),
                 ["Estimate","Std.Error","t value"],
                 ["x$i" for i = 1:length(cc)])
    end
end