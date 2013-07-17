# Author: Vincent Arel-Bundock
# Contact: varel@umich.edu
# Date: 2013-07-13
# License: BSD-3
# Original code from the python statsmodels project
# https://github.com/statsmodels/statsmodels

using DataFrames, Distributions

function bandwidth_hall_sheather(n, q, alpha)
    z = quantile(Normal(), q)
    num = 1.5 * pdf(Normal(), z)^2
    den = 2. * z^2 + 1
    h = n^(-1 / 3) * quantile(Normal(),(1. - alpha / 2.))^(2./3) * (num / den)^(1./3)
    h
end

function kernel_epanechnikov(x)
    u = copy(x)
    for i in 1:length(u)
        if abs(u[i]) <= 1
            u[i] = 3/4*(1-u[i]^2)
        else
            u[i] = 0
        end
    end
    u
end

function qreg_coef(y, X, q=.5)
    beta = ones(2)
    xstar = copy(X)
    diff = 10
    tol = 1e-6
    max_iter = 1000
    for i = 1:max_iter
        if diff > tol
            beta0 = copy(beta)
            xtx = *(xstar', X)
            xty = *(xstar', y)
            beta = *(pinv(xtx), xty)
            resid = y - *(X, beta)
            for j in 1:length(resid)
                if abs(resid[j]) < 1e-5
                    resid[j] = sign(resid[j]) * 1e-5
                end
                if resid[j] < 0
                    resid[j] = q * resid[j]
                else
                    resid[j] = (1-q) * resid[j]
                end
            end
            resid = abs(resid)
            xstar = broadcast(/, X, resid)
            #diff = max(abs(beta0-beta))
        end
    end
    beta[:,1]
end

function qreg_vcov(y, X, beta, q)
    resid = y - *(X, beta)
    iqre = quantile(resid[:,1], .75) - quantile(resid[:,1], .25)
    h = bandwidth_hall_sheather(length(resid), q, .05)
    h1 = min(std(y), iqre / 1.34) 
    h2 = quantile(Normal(), q+h) - quantile(Normal(), q-h)
    h = h1 * h2
    fhat0 = 1 / (length(resid)*h) * sum(kernel_epanechnikov(resid/h))
    d = copy(resid)
    for i in 1:length(d)
        if d[i] > 0
            d[i] = (q/fhat0)^2
        else
            d[i] = ((1-q)/fhat0)^2
        end
    end
    xtxi = pinv(*(X',X))
    xtdx = *(broadcast(*, X, d)', X)
    vcov = *(*(xtxi, xtdx), xtxi)
    vcov
end

function qreg(f::Expr, df::AbstractDataFrame, q::FloatingPoint=.5)
    mf = ModelFrame(f, df)
    mm = model_matrix(mf)
    mr = model_response(mf)
    coef = qreg_coef(mr.data, mm.m, q)
    vcov = qreg_vcov(mr.data, mm.m, coef, q)
    QRegModel(coef, vcov, mf)
end

type QRegModel
    beta::Array{Float64,1}
    vcov::Array{Float64}
    mf::ModelFrame
end

coef(x::QRegModel) = x.beta
vcov(x::QRegModel) = x.vcov
function summary(x::QRegModel)
    se = sqrt(diag(x.vcov))
    out = DataFrame()
    out["Estimate"] = x.beta
    out["Std.Error"] = se
    out["t value"] = x.beta ./ se
    out["2.5%"] = x.beta - 1.96 * se
    out["97.5%"] = x.beta + 1.96 * se
    out
end
