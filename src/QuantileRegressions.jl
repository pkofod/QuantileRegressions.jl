module QuantileRegressions

    import StatsModels: TableRegressionModel, FormulaTerm, coef, @formula
    export qreg, coef, vcov, stderr, quantiles, IP, IRLS, @formula

    using DataFrames, Distributions, LinearAlgebra, LinearAlgebra.BLAS, StatsModels, StatsBase

    mutable struct QRegModel
        beta::Vector{Float64}
        vcov::Matrix{Float64}
        stderr::Vector{Float64}
        q
    end

    abstract type Solver end

    struct IP <: Solver
    end

    struct IRLS <: Solver
        tol
        maxIter
        threshold
    end
    IRLS(;tol::Real = 1e-12, maxIter::Integer = 1_000, threshold::Real = 1e-5) = IRLS(tol, maxIter, threshold)

    include("InteriorPoint.jl")
    include("IRLS.jl")
    include("Covariance.jl")
    include("npqreg.jl")
    
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
        mm = ModelMatrix(mf)
        mr = response(mf)
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
