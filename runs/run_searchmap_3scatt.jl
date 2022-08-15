include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.)
normalized = false

## EP parameters ##
ϕinput = rand(3).+rand(3).*1im

ϵ = 0.1
## Evaluation ##
xx = LinRange(1/(4*sqrt(2)*π^2)-ϵ, 1/(4*sqrt(2)*π^2)+ϵ, 20) # gamma scan
yy = LinRange(0, 1000, 20) # 3rd scatt position scan
map = zeros((length(xx), length(yy)))
p = Progress(length(xx));
for i in eachindex(xx)
    Threads.@threads for j in eachindex(yy)
        scattpos = [[-π 0. 0.]
                    [0. yy[i] 0.] # changing the y location of the central scatterer does not change the resulting output power
                    [π 0. 0.]]
        alphas = [α(ω, 1, xx[i]) # it should be α(1, 1, 1/(4*sqrt(2)*π^2))
                  10^(4) # α(1, 1, 1000)
                  conj(α(ω, 1, xx[i]))]
        maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])
        map[j, i] = powerout(maxradius, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=0)[1]
    end
    next!(p)
end

npzwrite("./data/searchmap3EP.npz", Dict("xx" => xx, "yy" => yy, "map" => map))