include("./ScatterersDDA.jl")
using .ScatterersDDA
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.);

scattpos = [[0.5 0. 0.]
            [1. 0. 0.]
            [1.5 0. 0.]
            [2. 0. 0.]];
alphas = [(1. + 0im)
          (1. + 0im)
          (1. + 0im)
          (1. + 0im)];
ϕinput = [(1. + 0im)
          (1. + 0im)
          (1. + 0im)
          (1. + 0im)];

# scattpos = [0. 0. 0.];
# alphas = [(1. + 2im)];
# ϕinput = [(1. + 1im)];

rspan = LinRange(0, 4, 100)
Pout = zeros(length(rspan))
p = Progress(length(rspan));
Threads.@threads for i in 1:length(rspan)
    Pout[i] = powerout(rspan[i], ϕinput, scattpos, alphas, ω, J)[1]
    next!(p)
end

analyticalsum = zeros(length(scattpos[:, 1]))
for i in 1:length(scattpos[:, 1])
    analyticalsum[i] = evalsumm(i, ϕinput, scattpos, alphas, ω, J)
end

npzwrite("./data/data.npz", Dict("scattpos" => scattpos, "analyticalsum" => analyticalsum, "r" => rspan, "P" => Pout))