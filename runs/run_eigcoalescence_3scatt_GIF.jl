include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter

ω = 1.
J = Stdd(1.)
normalized = false

## EP parameters ##
ϵ = 0.02
rover2 = π
steps = 400
rspan = LinRange(sqrt(3)*rover2-ϵ, sqrt(3)*rover2+ϵ, steps)

###################
### Evaluated using Mathematica symbolic solution
### NB: there is a conj() in alphas[2] due to inconsistency with the matrix from mathematica.
###     With this conj() the 2 expressions coincide and we find the EP.
# G12 = 1/8/π^2
# alphas = [(1/2)*(1/G12 + ((3 + 2*sqrt(2))^(1/3)*Complex(G12^6)^(1/3))/G12^3 + G12/Complex(3*G12^6 + 2*sqrt(2)*G12^6)^(1/3))
#           conj(-((-2 + ((3 + 2*sqrt(2))^(1/3)*Complex(G12^6)^(1/3))/G12^2 + Complex(G12^6)^(2/3)/((3 + 2*sqrt(2))^(1/3)*G12^4)
#             + sqrt(6 - (3*Complex(G12^6)^(1/3))/((3 + 2*sqrt(2))^(2/3)*G12^2) - (3*Complex(3 + 2*sqrt(2))^(2/3)*Complex(G12^6)^(2/3))/G12^4))/(4*G12)))
#           (2 - ((3 + 2*sqrt(2))^(1/3)*Complex(G12^6)^(1/3))/G12^2 - G12^2/Complex(3*G12^6 + 2*sqrt(2)*G12^6)^(1/3) +
#             sqrt(6 - (3*(3 + 2*sqrt(2))^(2/3)*Complex(G12^6)^(2/3))/G12^4 - (3*G12^4)/Complex(3*G12^6 + 2*sqrt(2)*G12^6)^(2/3)))/(4*G12)]
###################
### Evaluated using NLsolve solution
alphas = [-7.013368475092600482445835849401992129215372596298039756222201423284825304069556 - 42.53024981045743992016289461947802809987344773660556638804706482651633259017354im, 132.4619897632574997546732604945833408650191728587364235821477266247433660309102 + 7.006613806219794923619144866639183660810213598907112437214055568579551121144819e-30im, -7.013368475092600482445835849402012558805883329658692400789161106750101841856979 + 42.53024981045743992016289461947723489761815542995479655135939519941658644818625im]
###################

## Evaluation ##
p = Progress(length(rspan));
eigs = zeros(Complex{Float64}, (steps, length(alphas)))
Threads.@threads for k in eachindex(rspan)
    scattpos = [[-rover2 0. 0.]
                [0. rspan[k] 0.] # equidistant scatterers with G = 1/8/π^2
                [rover2 0. 0.]]
    eigs[k, :] = eigvals(intmatrix(scattpos, alphas, ω, J; normalized=normalized, imagshift=1E-23))
    next!(p)
end

scattpos = [[-rover2 0. 0.]
            [0. rspan[1] 0.] # equidistant scatterers with G = 1/8/π^2
            [rover2 0. 0.]]
npzwrite("./data/gif_nlsolve.npz", Dict("rspan" => rspan, "rover2" => float(rover2), "epsilon" => ϵ, "alphas" => alphas, "scattpos" => scattpos, "eigs" => eigs))