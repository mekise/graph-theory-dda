include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.)
normalized = false

## spherical uniform scatt positions ##
# dim = 4
# maxradius = 5
# rscatt = rand(dim) .*maxradius
# cosθ = rand(dim) .*2 .-1
# sinθ = sqrt.(cosθ.^2)
# ϕ = rand(dim) .*2π

# scattx = rscatt .*cos.(ϕ) .*sinθ
# scatty = rscatt .*sin.(ϕ) .*sinθ
# scattz = rscatt .*cosθ
# scattpos = [scattx scatty scattz]

# alphas = (rand(dim) .+ rand(dim).*1im).*2 .-1 .-1im
# ϕinput = (rand(dim) .+ rand(dim).*1im).*2 .-1 .-1im

## in-line scatt positions ##
# scattpos = [[0.5 0. 0.]
#             [1. 0. 0.]
#             [1.5 0. 0.]
#             [2. 0. 0.]];
# alphas = [(1. + 0im)
#           (1. + 0im)
#           (1. + 0im)
#           (1. + 0im)];
# ϕinput = [(1. + 0im)
#           (1. + 0im)
#           (1. + 0im)
#           (1. + 0im)];
# maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

## single scatt tests ##
# scattpos = [0. 0. 0.];
# alphas = [(1. + 0im)];
# ϕinput = [(1. + 0im)];
# maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

## EP parameters ##
scattpos = [[-π 0. 0.]
            [0. 0. 0.] # changing the y location of the central scatterer does not change the resulting output power
            [π 0. 0.]]
alphas = [α(ω, 1, 1/(4*sqrt(2)*π^2)) # it should be α(1, 1, 1/(4*sqrt(2)*π^2))
          10^(-15) # α(1, 1, 1000)
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