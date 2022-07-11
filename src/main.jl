include("./ScatterersDDA.jl")
using .ScatterersDDA
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.);

scattpos = [[1. 2. 0.]
            [1. 1. 0.]
            [0. 2. 2.]];
alphas = [(1. + 2im)
          (1. - 2im)
          (1. - 1im)];
ϕinput = [(1. + 1im)
          (2. + 1im)
          (0. - 1im)];

rspan = LinRange(0, 4, 150)
Pout = zeros(length(rspan))

p = Progress(length(rspan));
Threads.@threads for i in 1:length(rspan)
    Pout[i] = powerout(rspan[i], ϕinput, scattpos, alphas, ω, J)[1]
    next!(p)
end

npzwrite("./data/power.npz", Dict("r" => rspan, "P" => Pout))
npzwrite("./data/scattpos.npz", Dict("scattpos" => scattpos))