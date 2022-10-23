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

r = 0.06
scattpos = zeros((nscatt, 3))
for j in 1:nscatt
    scattpos[j, :] = [r*cos(2π*j/nscatt) r*sin(2π*j/nscatt) 0.]
end

navg = 1000

# ϕinput = rand(nscatt).+rand(nscatt).*1im
ϕplane = r -> exp(1im*k0(ω, J)*r)
ϕinput  = zeros(ComplexF64, (nscatt))
for i in 1:nscatt
    global ϕinput[i] = ϕplane(scattpos[i, 1]) # plane wave propagating along the x-axis
end

# Total field
# xx = LinRange(-10, 10, 1000)
# phitot = zeros(ComplexF64, (length(xx)))
# weakphitot = zeros(ComplexF64, (length(xx)))
# strongphitot = zeros(ComplexF64, (length(xx)))
# for i in eachindex(xx)
#     phitot[i] = totfield([xx[i], 1., 0], ϕinput, scattpos, alphas, ω, J; imagshift=1E-23)# + ϕplane(xx[i])
#     weakphitot[i] = weaktotfield([xx[i], 1., 0], ϕinput, scattpos, alphas, ω, J; imagshift=1E-23)# + ϕplane(xx[i])
#     strongphitot[i] = strongtotfield([xx[i], 1., 0], ϕinput, scattpos, alphas, ω, J; imagshift=1E-23)# + ϕplane(xx[i])
# end
# npzwrite("./data/strong_approx.npz", Dict("r" => float(r), "alphas" => alphas, "scattpos" => scattpos, "phiinput" => ϕinput, "xx" => xx, "phitot" => phitot, "weakphitot" => weakphitot, "strongphitot" => strongphitot))

# Total field - map
xx = LinRange(-0.1, 0.1, 200)
yy = LinRange(-0.1, 0.1, 200)
phitot = zeros(ComplexF64, (length(xx), length(yy)))
strongphitot = zeros(ComplexF64, (length(xx), length(yy)))
deviationavg = zeros(Float64, (length(xx), length(yy)))
for k in 1:navg
    alphas = (rand(nscatt) .+ rand(nscatt).*1im)
    for i in eachindex(xx)
        for j in eachindex(yy)
            phitot[j, i] = totfield([xx[i], yy[j], 0], ϕinput, scattpos, alphas, ω, J; imagshift=1E-23)
            strongphitot[j, i] = strongtotfield([xx[i], yy[j], 0], ϕinput, scattpos, alphas, ω, J; imagshift=1E-23)
        end
    end
    deviationavg .+= abs.(phitot.-strongphitot)./abs.(phitot)
end
deviationavg .= deviationavg./navg

# npzwrite("./data/strong_approx_map.npz", Dict("r" => float(r), "alphas" => alphas, "scattpos" => scattpos, "phiinput" => ϕinput, "xx" => xx, "yy" => yy, "phitot" => phitot, "strongphitot" => strongphitot))
npzwrite("./data/strong_approx_avg.npz", Dict("r" => float(r), "scattpos" => scattpos, "phiinput" => ϕinput, "xx" => xx, "yy" => yy, "deviationavg" => deviationavg))