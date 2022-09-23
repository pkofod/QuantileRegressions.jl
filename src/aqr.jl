#  min τ*e'*u + (1-τ)*e'*v + λ*e'*w + λ*e'*z
#
#  s.t.
#       α - y  = v - u
#       D*α = z - w
#       v_i, u_i, z_i, w_i>=0
#       α_i free
#=
min sum_i(tau*u_i+(1-tau)*v_i) + sum_j(x_j+z_j)

y-alpha_i = u_i - v_i
d_j'alpha = x_j-z_j, j = 1...n-1
u_i >= 0
v_i >= 0
x_j >= 0
z_j >= 0
=#

using LinearAlgebra
using SparseArrays
using QuadraticModels
using Plots
using Distributions
using RipQP
using Random
using DataInterpolations

struct AQR
end

function train_eval_split(X, Y, n, rng=Random.GLOBAL_RNG)
    excl = rand(rng, 1:length(Y), n)
    incl = [i for i in 1:length(Y) if i ∉ excl]
    Xtrain = X[incl]
    Ytrain = Y[incl]
    Xeval = X[excl]
    Yeval = Y[excl]
    Xtrain, Ytrain, Xeval, Yeval
end

function nonuniqueI(xsorted)
	_x = first(xsorted)
	VALS = eltype(xsorted)[]
	COLS = Int[]
	ROWS = Int[]
	j = 1
	i = 1
for __x in xsorted
    if _x != __x
        _x = __x
        j += 1
    end
    push!(VALS, 1)
    push!(COLS, j)
    push!(ROWS, i)
    i += 1
end
sparse(ROWS, COLS, VALS)
end

function rqss(x, y, τ, λ, evalx=nothing, evaly=nothing)
    xperm = sortperm(x)
    xsorted = x[xperm]
    xknots = unique(xsorted)
    ysorted = y[xperm]
    xdelta = diff(xknots)
    xdeltainv = inv.(xdelta)
    N = length(x)
    n = length(xknots)
    D = zeros(n, n)
    for j=1:n-2
        D[j,j] = xdeltainv[j]
        D[j,j+1] = -xdeltainv[j+1]-xdeltainv[j]
        D[j,j+2] = xdeltainv[j+1]
    end
    D
    DI = Matrix(I, size(D)...)
    D0 = spzeros(N, n)
    XI = nonuniqueI(xsorted)
    XI0 = spzeros(n, N)
    XIy = sparse(I, N, N)

    α = zeros(n)
    Design1 = hcat(XI, XIy, -XIy, D0,  D0)
    Design2 = hcat(D,  XI0,  XI0, DI, -DI)
    A = vcat(Design1, Design2)
    c = vcat(
        zeros(n), # coefficients
        τ*ones(N), # u
        (1-τ)*ones(N), # v
        λ*ones(n), # penalty w
        λ*ones(n), # penalty z
        )
    rhs = vcat(ysorted, zeros(n))
    l=vcat(fill(-Inf, n),
           zeros(N),
           zeros(N),
           zeros(n),
           zeros(n) # coef
            )

    QM = QuadraticModel(c, spzeros(length(c),length(c)), A=A, lcon=rhs, ucon=rhs, lvar=l,c0=0.0, name="RQSS")
    sol = ripqp(QM; display=false)

    interp = DataInterpolations.LinearInterpolation(sol.solution[1:length(xknots)], xknots)
    if !(evalx === nothing && evaly === nothing)
        evalresid = evaly.-interp.(evalx)
    else
        evalresid = nothing
    end
    (; τ, λ, evalresid, f=interp)
end


function aqr(X, Y, τ, lambda)
	#=
	LAM = []
	Λ = [10.0^i for i = -9:3]
	τ = 0.5
	for l = Λ
	     push!(LAM, rqss(X, Y, τ, l))
	end
	plt=Plots.scatter(X, Y, label="data", mc=:black, ms=2)
	for i = eachindex(LAM)
	    Plots.plot!(sort(unique(X)), LAM[i].f.(sort(unique(X))); label="τ=$(τ), λ=$(LAM[i].λ)", lw=2, alpha=0.9)
	end
	plt
	=#
	if lambda === nothing
		Λ = vcat([0.001,0.01, 0.05,0.1,0.2, 0.25,0.3,0.5, 0.75], [1.0, 2.0,5.0])
		absmaxdiff = abs(-(extrema(Y)...))*100
		for i = 1:50
			teni = 10.0^i
			if teni < absmaxdiff
				push!(Λ, teni)
			end
		end
		sse = zeros(length(Λ))
		for i = 1:80
		    Xtrain, Ytrain, Xeval, Yeval = train_eval_split(X, Y, floor(Int, length(Y)*0.8))
		    for j = 1:length(Λ)
		        for k in [0.5]#[0.1,0.5,0.9]
		           res = rqss(Xtrain, Ytrain, k, Λ[j], Xeval, Yeval)
		           sse[j] += sum(x->x^2, res.evalresid)
		        end
		    end
		end
		best = rqss(X, Y, τ, Λ[findmin(sse)[2]])
	else
		rqss(X, Y, τ, lambda)
	end
#	best = rqss(X, Y, τ, 2.0)
end
#=
Plots.plot(Λ, sqrt.(sse))
findmin(sse)
function plot_best(sse, X, Y, τ)
    plt=Plots.scatter(X, Y, label="data", mc=:black, ms=2)
    for t ∈ τ
        best = rqss(X, Y, t, Λ[findmin(sse)[2]])
        Plots.plot!(sort(unique(X)), best.f.(sort(unique(X))); label="τ=$(t), λ=$(Λ[findmin(sse)[2]])", lw=2, alpha=0.9)
    end
    plt
end
function plot_specific(sse, X, Y, τ, i)
    plt=Plots.scatter(X, Y, label="data", mc=:black, ms=2)
    for t ∈ τ
        best = rqss(X, Y, t, Λ[i])
        Plots.plot!(sort(unique(X)), best.f.(sort(unique(X))); label="τ=$(t), λ=$(Λ[i])", lw=2, alpha=0.9)
    end
    plt
end
plot_best(sse, X, Y, [0.1, 0.5, 0.9])
=#