include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.)
normalized = false

## EP parameters ##
ϵ = 0.05
ωspan = LinRange(ω-ϵ, ω+ϵ, 100)
scattpos = [[-π 0. 0.]
            [0. 0. 0.] # changing the y location of the central scatterer does not change the resulting output power
            [π 0. 0.]]
ϕinput = rand(3).+rand(3).*1im
maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

## Evaluation ##
xx = LinRange(-10, 10, 200)
yy = LinRange(-10, 10, 200)
phitot = zeros(ComplexF64, (length(xx), length(yy)))
p = Progress(length(xx));
for k in eachindex(ωspan)
    for i in eachindex(xx)
        Threads.@threads for j in eachindex(yy)
            alphas = [α(ωspan[k], 1, 1/(4*sqrt(2)*π^2)) # it should be α(1, 1, 1/(4*sqrt(2)*π^2))
                      α(ωspan[k], 1, 0.00000001)
                      conj(α(ωspan[k], 1, 1/(4*sqrt(2)*π^2)))]
            phitot[j, i] = totfield([xx[i], yy[j], 0], ϕinput, scattpos, alphas, ω, J; normalized=normalized)
        end
    end
    next!(p)
    npzwrite("./gif/gif"*string(k)*".npz", Dict("scattpos" => scattpos, "xx" => xx, "yy" => yy, "phitot" => phitot))
end