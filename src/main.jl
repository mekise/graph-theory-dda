include("./ScatterersDDA.jl")
using .ScatterersDDA
using NPZ

ω = 1.
J = Stdd(1.);

scattpos = [[0. 0. 0.]
            [1. 1. 1.]];
alphas = [(1. + 2im)
          (1. - 2im)];
ϕinput = [(1. +1im)
          (2. +1im)];

# scattpos = [[1. 2. 3.]];
# alphas = [(1. + 2im)];
# ϕinput = [(1. +1im)];

# rspan = LinRange(0, 10, 10)
# Pout = zeros(length(rspan))
# Threads.@threads for i in 1:length(rspan)
#     Pout[i] = powerout(rspan[i], ϕinput, scattpos, alphas, ω, J)[1]
# end

# npzwrite("./data/power.npz", Dict("r" => rspan, "P" => Pout))
# npzwrite("./data/scattpos.npz", Dict("scattpos" => scattpos))

res = Gintegral(2, ω, J)

npzwrite("./data/Gintegral.npz", Dict("res" => res))