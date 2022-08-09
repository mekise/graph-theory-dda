include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

J = Stdd(1.)
normalized = false

## EP parameters ##
scattpos = [[-π 0. 0.]
            [0. 0. 0.] # changing the y location of the central scatterer does not change the resulting output power
            [π 0. 0.]]
ϕinput = rand(3).+rand(3).*1im
maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

## Evaluation ##
ϵ = 0.1
ωspan = LinRange(1-ϵ, 1+ϵ, 100)
Pout = zeros(length(ωspan))
Poutexpl = zeros(length(ωspan))
p = Progress(length(ωspan))
println("Starting...")
Threads.@threads for i in eachindex(ωspan)
    alphas = [α(ωspan[i], 1, 1/(4*sqrt(2)*π^2)) # it should be α(1, 1, 1/(4*sqrt(2)*π^2))
              10^(-15) # α(1, 1, 1000)
              conj(α(ωspan[i], 1, 1/(4*sqrt(2)*π^2)))]
    Pout[i] = powerout(maxradius, ϕinput, scattpos, alphas, ωspan[i], J; normalized=normalized, imagshift=0)[1]
    # Poutexpl[i] = poweroutexplicit(maxradius, ϕinput, scattpos, alphas, ωspan[i], J; normalized=normalized)[1]
    next!(p)
end

npzwrite("./data/poweroveromega3EP.npz", Dict("scattpos" => scattpos, "omega" => ωspan, "P" => Pout))