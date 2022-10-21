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

alphas = (rand(nscatt) .+ rand(nscatt).*1im)

r = 2.
scattpos = zeros((nscatt, 3))
for j in 1:nscatt
    scattpos[j, :] = [r*cos(2π*j/nscatt) r*sin(2π*j/nscatt) 0.]
end

# ϕinput = rand(nscatt).+rand(nscatt).*1im
ϕplane = r -> exp(1im*k0(ω, J)*r)
ϕinput  = zeros(ComplexF64, (nscatt))
for i in 1:nscatt
    global ϕinput[i] = ϕplane(scattpos[i, 1]) # plane wave propagating along the x-axis
end

# Total field
xx = LinRange(-10, 10, 1000)
phitot = zeros(ComplexF64, (length(xx)))
for i in eachindex(xx)
    phitot[i] = totfield([xx[i], 0.5, 0], ϕinput, scattpos, alphas, ω, J; imagshift=1E-23) + ϕplane(xx[i])
end
npzwrite("./data/approximation.npz", Dict("r" => float(r), "alphas" => alphas, "scattpos" => scattpos, "phiinput" => ϕinput, "xx" => xx, "phitot" => phitot))