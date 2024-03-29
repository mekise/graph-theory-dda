include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

J = Stdd(1.)
normalized = true

## EP parameters ##
scattpos = [[0. 0. 0.]
            [2π 0. 0.]]
ϕinput = rand(2).+rand(2).*1im
maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

## Evaluation ##
ϵ = 0.1
ωspan = LinRange(1-ϵ, 1+ϵ, 200)
Pout = zeros(length(ωspan))
Poutexpl = zeros(length(ωspan))
p = Progress(length(ωspan));
Threads.@threads for i in eachindex(ωspan)
    alphas = [α(ωspan[i], 1, 1/(8*π^2))
              -α(ωspan[i], 1, 1/(8*π^2))]
    Pout[i] = powerout(maxradius, ϕinput, scattpos, alphas, ωspan[i], J; normalized=normalized)[1]
    Poutexpl[i] = poweroutexplicit(maxradius, ϕinput, scattpos, alphas, ωspan[i], J; normalized=normalized)[1]
    next!(p)
end

npzwrite("./data/poweroveromeganorm.npz", Dict("scattpos" => scattpos, "omega" => ωspan, "P" => Pout, "Pexpl" => Poutexpl))