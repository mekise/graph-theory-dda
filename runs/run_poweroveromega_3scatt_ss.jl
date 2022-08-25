include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

J = Stdd(1.)
normalized = false

# this is the set of parameters to reproduce my matrix for higher order EPs

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
ϵ = 0.1
ωspan = LinRange(1-ϵ, 1+ϵ, 1000)
Pout = zeros(length(ωspan))
p = Progress(length(ωspan))

println("Starting...")
Threads.@threads for i in eachindex(ωspan)
    Pout[i] = powerout(maxradius, ϕinput, scattpos, alphas, ωspan[i], J; normalized=normalized)[1]
    next!(p)
end

npzwrite("./data/poweroveromega3EPss.npz", Dict("scattpos" => scattpos, "omega" => ωspan, "P" => Pout))