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

scattpos = zeros((steps, 3, 3))
for i in 1:steps
    scattpos[i, :, :] = [[-rover2 0. 0.]
                         [0. rspan[i] 0.] # equidistant scatterers with G = 1/8/π^2
                         [rover2 0. 0.]]
end
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

alphas = [132.46198976325752 + 0.0im
-7.0133684750926255 - 42.53024981045747im
  -7.01336847509261 + 42.53024981045747im]
###################
### Evaluated using NLsolve solution
# alphas = [-7.013368475092600482445835849401586148371126697377497396993147566773778372522662 + 42.53024981045743992016289461947781909647920721596898054239751469515387324972564im, -7.013368475092600482445835849401586148371126697377497396993147566773778372522731 - 42.53024981045743992016289461947781909647920721596898054239751469515387324972564im, 132.46198976325749975467326049460771358600803016240562244219341392091438694475 + 6.703134469499830002807082037260090296948103626424876531613289800474263330871294e-77im]
###################

## Evaluation ##
p = Progress(length(rspan));
eigs = zeros(Complex{Float64}, (steps, length(alphas)))
Threads.@threads for k in eachindex(rspan)
    eigs[k, :] = eigvals(intmatrix(scattpos[k, :, :], alphas, ω, J; normalized=normalized, imagshift=1E-23))
    next!(p)
end

npzwrite("./data/eigcoalescence_3scatt_test.npz", Dict("rspan" => rspan, "rover2" => float(rover2), "epsilon" => ϵ, "alphas" => alphas, "scattpos" => scattpos, "eigs" => eigs))