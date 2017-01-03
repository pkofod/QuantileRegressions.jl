using DataFrames

module QuantileRegression

    import DataFrames.DataFrameRegressionModel
    import DataFrames.coef
    import Base.LinAlg.BLAS.axpy!
    export qreg, coef, vcov, stderr, quantiles, IP, IRLS

    using DataFrames, Distributions, Base.LinAlg.BLAS


    type QRegModel
        beta::Vector{Float64}
        vcov::Matrix{Float64}
        stderr::Vector{Float64}
        q
    end

    abstract Solver

    immutable IP <: Solver
    end

    immutable IRLS <: Solver
    end

    include("InteriorPoint.jl")
    include("IRLS.jl")
    include("Covariance.jl")

    function qreg(f::Formula, df::AbstractDataFrame, q, s::Solver = IP())
        mf = ModelFrame(f, df)
        mm = ModelMatrix(mf)
        mr = model_response(mf)
        coef = qreg_coef(mr, mm.m, q, s)
        vcov = qreg_vcov(mr, mm.m, coef, q)
        stderr = sqrt(diag(vcov))
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
