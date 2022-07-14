include("./ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.);

## spherical uniform scatt positions ##
dim = 4
maxradius = 5
rscatt = rand(dim) .*maxradius
cosθ = rand(dim) .*2 .-1
sinθ = sqrt.(cosθ.^2)
ϕ = rand(dim) .*2π
scattx = rscatt .*cos.(ϕ) .*sinθ
scatty = rscatt .*sin.(ϕ) .*sinθ
scattz = rscatt .*cosθ
scattpos = [scattx scatty scattz]
alphas = (rand(dim) .+ rand(dim).*1im).*2 .-1 .-1im
ϕinput = (rand(dim) .+ rand(dim).*1im).*2 .-1 .-1im

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

## single scatt tests ##
# scattpos = [0. 0. 0.];
# alphas = [(1. + 2im)];
# ϕinput = [(1. + 1im)];

# maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])
rspan = LinRange(0, maxradius + 1, 50)
Pout = zeros(length(rspan))
p = Progress(length(rspan));
Threads.@threads for i in 1:length(rspan)
    Pout[i] = powerout(rspan[i], ϕinput, scattpos, alphas, ω, J; imagshift=0)[1]
    next!(p)
end

analyticalsum = zeros(length(scattpos[:, 1]))
for i in 1:length(scattpos[:, 1])
    analyticalsum[i] = evalsumm(i, ϕinput, scattpos, alphas, ω, J)
end

npzwrite("./data/data.npz", Dict("scattpos" => scattpos, "analyticalsum" => analyticalsum, "r" => rspan, "P" => Pout))