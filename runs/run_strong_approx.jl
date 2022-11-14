include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter
using NLsolve
using Combinatorics

ω = 1.
J = Stdd(1.)
nscatt = 4

# alphas = (rand(nscatt) .+ rand(nscatt).*1im)

r = 0.001
scattpos = zeros((nscatt, 3))
for j in 1:nscatt
    scattpos[j, :] = [r*cos(2π*j/nscatt) r*sin(2π*j/nscatt) 0.]
end

navg = 100

# ϕinput = rand(nscatt).+rand(nscatt).*1im
ϕplane = r -> exp(1im*k0(ω, J)*r)
ϕinput  = zeros(ComplexF64, (nscatt))
for i in 1:nscatt
    global ϕinput[i] = ϕplane(scattpos[i, 1]) # plane wave propagating along the x-axis
end

# Total field
xx = LinRange(-1.5, 1.5, 200)
yy = LinRange(-1.5, 1.5, 200)
phitot = zeros(ComplexF64, (length(xx), length(yy)))
phitotapprox = zeros(ComplexF64, (length(xx), length(yy)))
deviationavg = zeros(Float64, (length(xx), length(yy)))
for k in 1:navg
    alphas = (rand(nscatt) .+ rand(nscatt).*1im)
    for i in eachindex(xx)
        for j in eachindex(yy)
            phitot[j, i] = totfield([xx[i], yy[j], 0], ϕinput, scattpos, alphas, ω, J; imagshift=1E-23)
            phitotapprox[j, i] = strongtotfield([xx[i], yy[j], 0], ϕinput, scattpos, alphas, ω, J; imagshift=1E-23)
        end
    end
    deviationavg .+= abs.(phitot.-phitotapprox) ./ abs.(phitot)
end
deviationavg ./= navg

npzwrite("./data/strong_approx_map.npz", Dict("r" => float(r), "scattpos" => scattpos, "phiinput" => ϕinput, "xx" => xx, "yy" => yy, "deviationavg" => deviationavg, "phitot" => phitot, "strongphitot" => phitotapprox))