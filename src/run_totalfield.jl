include("./ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.);

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
scattpos = [[0. 0. 0.]
            [2π 0. 0.]]
alphas = [α(1, 1, 1/(8*π^2))
          -α(1, 1, 1/(8*π^2))]
ϕinput = rand(2).+rand(2).*1im
maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

## Evaluation ##
xx = LinRange(-10, 10, 200)
yy = LinRange(-10, 10, 200)
phitot = zeros(ComplexF64, (length(xx), length(yy)))
p = Progress(length(xx));
for i in 1:length(xx)
    Threads.@threads for j in 1:length(yy)
        phitot[j, i] = totfield([xx[i], yy[j], 0], ϕinput, scattpos, alphas, ω, J)
    end
    next!(p)
end

npzwrite("./data/data_total_field.npz", Dict("scattpos" => scattpos, "xx" => xx, "yy" => yy, "phitot" => phitot))