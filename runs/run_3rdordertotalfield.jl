include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.)
normalized = false

## EP parameters ##
scattpos = [[-π 0. 0.]
            [0. 0. 0.] # changing the y location of the central scatterer does not change the resulting output power
            [π 0. 0.]]
alphas = [α(ω, 1, 1/(4*sqrt(2)*π^2)) # it should be α(1, 1, 1/(4*sqrt(2)*π^2))
          α(ω, 1, 0.00000001)
          conj(α(ω, 1, 1/(4*sqrt(2)*π^2)))]
ϕinput = rand(3).+rand(3).*1im
maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

## Evaluation ##
xx = LinRange(-10, 10, 200)
yy = LinRange(-10, 10, 200)
phitot = zeros(ComplexF64, (length(xx), length(yy)))
p = Progress(length(xx));
for i in eachindex(xx)
    Threads.@threads for j in eachindex(yy)
        phitot[j, i] = totfield([xx[i], yy[j], 0], ϕinput, scattpos, alphas, ω, J; normalized=normalized)
    end
    next!(p)
end

npzwrite("./data/totalfield3EP.npz", Dict("scattpos" => scattpos, "xx" => xx, "yy" => yy, "phitot" => phitot))