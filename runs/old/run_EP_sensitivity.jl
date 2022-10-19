include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using ForwardDiff
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.)
normalized = false

function sensitivityEP(r, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-23, reltol=1e-8)
    f = ω -> powerout(r, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift, reltol=reltol)
    der = ω -> ForwardDiff.derivative(f, ω)
    return der(ω)
end

## EP parameters ##
scattpos = [[0. 0. 0.]
            [2π 0. 0.]]
alphas = [α(1, 1, 1/(8*π^2))
          -α(1, 1, 1/(8*π^2))]
ϕinput = rand(2).+rand(2).*1im
maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

## Evaluation ##
rspan = LinRange(0, maxradius*1.4, 100)
Pout = zeros(length(rspan))
Poutexpl = zeros(length(rspan))
p = Progress(length(rspan));
Threads.@threads for i in eachindex(rspan)
    Pout[i] = powerout(rspan[i], ϕinput, scattpos, alphas, ω, J; normalized=normalized)[1]
    Poutexpl[i] = poweroutexplicit(rspan[i], ϕinput, scattpos, alphas, ω, J; normalized=normalized)[1]
    next!(p)
end

analyticalsum = zeros(length(scattpos[:, 1]))
analyticalsumcorrected = zeros(length(scattpos[:, 1]))
for i in 1:length(scattpos[:, 1])
    analyticalsum[i] = evalsumm(i, ϕinput, scattpos, alphas, ω, J; normalized=normalized)
    analyticalsumcorrected[i] = evalsummcorrected(i, ϕinput, scattpos, alphas, ω, J; normalized=normalized)
end

npzwrite("./data/sensitivityEP.npz", Dict("scattpos" => scattpos, "analyticalsum" => analyticalsum, "analyticalsumcorrected" => analyticalsumcorrected, "r" => rspan, "P" => Pout, "Pexpl" => Poutexpl))