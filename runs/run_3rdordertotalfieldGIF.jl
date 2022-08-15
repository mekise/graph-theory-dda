include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.)
normalized = false

## EP parameters ##
ϵ = 0.05
ωspan = LinRange(ω-ϵ, ω+ϵ, 200)

## EP parameters ##
scattpos = [[-π 0. 0.]
            [0. sqrt(3)*π 0.] # equidistant scatterers with G = 1/8/π^2
            [π 0. 0.]]
ϕinput = rand(3).+rand(3).*1im
maxradius = maximum([norm(scattpos[i, :]) for i in 1:length(scattpos[:, 1])])

###################
### NB: there is a conj() in alphas[2] due to inconsistency with the matrix from mathematica.
###     With this conj() the 2 expressions coincide and we find the EP.
G12 = 1/8/π^2
alphas = [(1/2)*(1/G12 + ((3 + 2*sqrt(2))^(1/3)*Complex(G12^6)^(1/3))/G12^3 + G12/Complex(3*G12^6 + 2*sqrt(2)*G12^6)^(1/3))
          conj(-((-2 + ((3 + 2*sqrt(2))^(1/3)*Complex(G12^6)^(1/3))/G12^2 + Complex(G12^6)^(2/3)/((3 + 2*sqrt(2))^(1/3)*G12^4)
            + sqrt(6 - (3*Complex(G12^6)^(1/3))/((3 + 2*sqrt(2))^(2/3)*G12^2) - (3*Complex(3 + 2*sqrt(2))^(2/3)*Complex(G12^6)^(2/3))/G12^4))/(4*G12)))
          (2 - ((3 + 2*sqrt(2))^(1/3)*Complex(G12^6)^(1/3))/G12^2 - G12^2/Complex(3*G12^6 + 2*sqrt(2)*G12^6)^(1/3) +
            sqrt(6 - (3*(3 + 2*sqrt(2))^(2/3)*Complex(G12^6)^(2/3))/G12^4 - (3*G12^4)/Complex(3*G12^6 + 2*sqrt(2)*G12^6)^(2/3)))/(4*G12)]
###################

## Evaluation ##
xx = LinRange(-10, 10, 200)
yy = LinRange(-10, 10, 200)
phitot = zeros(ComplexF64, (length(xx), length(yy)))
p = Progress(length(xx));
for k in eachindex(ωspan)
    for i in eachindex(xx)
        Threads.@threads for j in eachindex(yy)
            phitot[j, i] = totfield([xx[i], yy[j], 0], ϕinput, scattpos, alphas, ωspan[k], J; normalized=normalized)
        end
    end
    next!(p)
    npzwrite("./gif/gif"*string(k)*".npz", Dict("omegaspan" => ωspan, "phiinput" => ϕinput, "scattpos" => scattpos, "xx" => xx, "yy" => yy, "phitot" => phitot))
end