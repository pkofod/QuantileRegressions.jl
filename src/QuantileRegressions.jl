module QuantileRegressions

    import StatsModels: TableRegressionModel, FormulaTerm, coef, @formula
    export qreg, coef, vcov, stderr, quantiles, IP, IRLS, @formula
    import LinearAlgebra.BLAS
    using DataFrames, Distributions, LinearAlgebra, StatsModels, StatsBase

    mutable struct QRegModel
        beta::Vector{Float64}
        vcov::Matrix{Float64}
        stderr::Vector{Float64}
        q
    end

    abstract type Solver end

    struct IRLS{T} <: Solver
        tol::T
        threshold::T
        maxiter::Integer
    end
    IRLS(;tol::Real = 1e-12, maxiter::Integer = 1_000, threshold::Real = 1e-5) = IRLS(tol, threshold, maxiter)

    include("InteriorPoint.jl")
    include("IRLS.jl")
    include("Covariance.jl")
    include("npqreg.jl")
    include("aqr.jl")

    function qreg(f::FormulaTerm, df::AbstractDataFrame, q, weights::AbstractVector, s::Solver = IP())
        mf = ModelFrame(f, df)
        mm = ModelMatrix(mf)
        mm = ModelMatrix(mm.m.*weights, mm.assign)
        mr = response(mf).*weights

        coef = qreg_coef(mr, mm.m, q, s)
        vcov = qreg_vcov(mr, mm.m, coef, q)
        stderror = sqrt.(diag(vcov))
        return TableRegressionModel(QRegModel(coef, vcov, stderror, q), mf, mm)
    end
    function qreg(f::FormulaTerm, df::AbstractDataFrame, q, s::Solver = IP())
        mf = ModelFrame(f, df)
        mm = ModelMatrix(mf)::ModelMatrix{Matrix{Float64}}
        mr = response(mf)::Vector{Float64}
        coef = qreg_coef(mr, mm.m, q, s)
        vcov = qreg_vcov(mr, mm.m, coef, q)
        stderror = sqrt.(diag(vcov))
        return TableRegressionModel(QRegModel(coef, vcov, stderror, q), mf, mm)
    end
    qreg(f::FormulaTerm, df::AbstractDataFrame, s::Solver) = qreg(f, df, 0.5, s)

    StatsBase.coef(x::QRegModel) = x.beta
    StatsBase.coef(rm::TableRegressionModel) = coef(rm.model)

    StatsBase.vcov(x::QRegModel) = x.vcov

    StatsBase.stderror(x::QRegModel) = x.stderr

    quantiles(x::QRegModel) = x.q
    quantiles(rm::TableRegressionModel) = quantiles(rm.model)

    function StatsBase.coeftable(mm::QRegModel)
        cc = coef(mm)
        se = stderror(mm)
        tt = cc./se
        CoefTable(hcat(kron(mm.q,ones(length(se))), cc,se,tt),
                 ["Quantile",  "Estimate","Std.Error","t value"],
                 ["x$i" for i = 1:length(cc)])
    end
end
