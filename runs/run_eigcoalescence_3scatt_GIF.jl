include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter
using NLsolve

ω = 1.
J = Stdd(1.)
normalized = false

###################
### EP parameters
ϵ = 0.02
rover2 = π
steps = 400
rspan = LinRange(sqrt(3)*rover2-ϵ, sqrt(3)*rover2+ϵ, steps)

scattpos = zeros((steps, 3, 3))
for i in 1:steps
    scattpos[i, :, :] = [[-rover2 0. 0.]
                         [0. rspan[i] 0.]
                         [rover2 0. 0.]]
end

###################
### Mathematica solutions 
### full
# alphas = [132.46198976325752 + 0.0im
#           -7.0133684750926255 - 42.53024981045747im
#           -7.01336847509261 + 42.53024981045747im]
### truncated
# alphas = [132.46 + 0.0im
#           -7.01 - 42.53im
#           -7.01 + 42.53im]

###################
### non-linear system solver
G12 = 1/8/π^2

function f!(F, x)
    F[1] = -G12^3 + 1/(2*x[1]*x[2]*x[3])
    F[2] = 3*G12^2 - (x[1]+x[2]+x[3])/(x[1]*x[2]*x[3])
    F[3] = x[1] + x[2]*x[3]/(x[2]+x[3])
end
function j!(J, x)
    J[1, 1] = -1/(2*x[1]^2*x[2]*x[3])
    J[1, 2] = -1/(2*x[1]*x[2]^2*x[3])
    J[1, 3] = -1/(2*x[1]*x[2]*x[3]^2)
    J[2, 1] = (x[2]+x[3])/(x[1]^2*x[2]*x[3])
    J[2, 2] = (x[1]+x[3])/(x[1]*x[2]^2*x[3])
    J[2, 3] = (x[1]+x[2])/(x[1]*x[2]*x[3]^2)
    J[3, 1] = 1
    J[3, 2] = x[3]/(x[2]+x[3]) - x[2]*x[3]/(x[2]+x[3])^2
    J[3, 3] = x[2]/(x[2]+x[3]) - x[2]*x[3]/(x[2]+x[3])^2
end
validsolution = "n"
while validsolution != "y"
    global s = nlsolve(f!, j!, 10*(rand(3)+rand(3)*1im), ftol=1e-20, iterations=1_000)
    while s.f_converged == false
        s = nlsolve(f!, j!, 10*(rand(3)+rand(3)*1im), ftol=1e-20, iterations=1_000)
    end
    s = nlsolve(f!, j!, s.zero, ftol=1e-80, iterations=10_000)
    FF = [1. + 1im, 1., 1.]
    f!(FF, s.zero)
    println("sol = ", s.zero)
    println("\nftol = " * string(s.ftol) * "\nresidual_norm = " * string(s.residual_norm) * "\nFF(alphas) = " * string(FF))
    println("\nContinue with this solution? (y, n)")
    global validsolution = readline()
end
alphas = s.zero

###################
### Evaluation
p = Progress(length(rspan));
eigs = zeros(Complex{Float64}, (steps, length(alphas)))
Threads.@threads for k in eachindex(rspan)
    eigs[k, :] = eigvals(intmatrix(scattpos[k, :, :], alphas, ω, J; normalized=normalized, imagshift=1E-23))
    next!(p)
end

npzwrite("./data/eigcoalescence_3scatt_nlsolve.npz", Dict("rspan" => rspan, "rover2" => float(rover2), "epsilon" => ϵ, "alphas" => alphas, "scattpos" => scattpos, "eigs" => eigs))