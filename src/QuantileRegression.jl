using DataFrames

module QuantileRegression

    import  DataFrames.DataFrameRegressionModel
    export qreg, coef, vcov, stderr

    using DataFrames, Distributions, Base.LinAlg.BLAS

    include("InteriorPoint.jl")
    include("IRLS.jl")
    include("Covariance.jl")

    type QRegModel
        beta::Vector{Float64}
        vcov::Matrix{Float64}
        stderr::Vector{Float64}
    end

    function qreg(f::Formula, df::AbstractDataFrame; q::Real = 0.5, method::Symbol = :ip)
        mf = ModelFrame(f, df)
        mm = ModelMatrix(mf)
        mr = model_response(mf)
        if method == :irls
            coef = qreg_coef(mr, mm.m, q)
        elseif method == :ip
            coef = qreg_ip_coef(mr, mm.m,q)
        else
            println("You didn't provide a supported method, defaulting to :ip.")
            coef = qreg_ip_coef(mr, mm.m, q)
        end
        vcov = qreg_vcov(mr, mm.m, coef, q)
        stderr = sqrt(diag(vcov))
        return DataFrameRegressionModel(QRegModel(coef, vcov, stderr),mf,mm)
    end

    coef(x::QRegModel) = x.beta

    vcov(x::QRegModel) = x.vcov

    stderr(x::QRegModel) = x.stderr

    function StatsBase.coeftable (mm::QRegModel)
        cc = coef(mm)
        se = stderr(mm)
        tt = cc./se
        CoefTable(hcat(cc,se,tt),
                 ["Estimate","Std.Error","t value"],
                 ["x$i" for i = 1:length(cc)])
    end
end
