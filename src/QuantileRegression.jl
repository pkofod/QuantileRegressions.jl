using DataFrames

module QuantileRegression

    import StatsModels: DataFrameRegressionModel, Formula, coef, @formula
    import LinearAlgebra: BLAS.axpy!
    import Statistics.quantile
    export qreg, coef, vcov, stderr, quantiles, IP, IRLS, @formula

    using DataFrames, Distributions, LinearAlgebra, StatsModels, StatsBase, Statistics

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

    function qreg(f::Formula, df::AbstractDataFrame, q, s::Solver = IP())
        mf = ModelFrame(f, df)
        mm = ModelMatrix(mf)
        mr = model_response(mf)
        coef = qreg_coef(mr, mm.m, q, s)
        vcov = qreg_vcov(mr, mm.m, coef, q)
        stderr = sqrt.(diag(vcov))
        return DataFrameRegressionModel(QRegModel(coef, vcov, stderr, q), mf, mm)
    end
    qreg(f::Formula, df::AbstractDataFrame, s::Solver) = qreg(f, df, 0.5, s)

    coef(x::QRegModel) = x.beta
    coef(rm::DataFrameRegressionModel) = coef(rm.model)

    vcov(x::QRegModel) = x.vcov

    stderr(x::QRegModel) = x.stderr

    quantiles(x::QRegModel) = x.q
    quantiles(rm::DataFrameRegressionModel) = quantiles(rm.model)

    function StatsBase.coeftable(mm::QRegModel)
        cc = coef(mm)
        se = stderr(mm)
        tt = cc./se
        CoefTable(hcat(kron(mm.q,ones(length(se))), cc,se,tt),
                 ["Quantile",  "Estimate","Std.Error","t value"],
                 ["x$i" for i = 1:length(cc)])
    end
end
