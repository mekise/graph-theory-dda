include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter
using NLsolve

ω = 1.
J = Stdd(1.)
normalized = false
nscatt = 4
steps = 400

###################
### EP parameters
ϵ = 0.02
r = π
diagover2 = sqrt(2)*r/2
rspan = LinRange(diagover2-ϵ, diagover2+ϵ, steps)

scattpos = zeros((steps, nscatt, 3))
for i in 1:steps
    scattpos[i, :, :] = [[diagover2 0. 0.]
                         [0. rspan[i] 0.]
                         [-diagover2 0. 0.]
                         [0. -diagover2 0.]]
end

###################
### non-linear system solver
G1 = greensfun(scattpos[1, 1, :], scattpos[1, 4, :], ω, J)
G2 = greensfun(scattpos[1, 1, :], scattpos[1, 3, :], ω, J)
function f!(F, x)
    F[1] = -G1^2*(x[1]+x[3]+2*G2*x[1]*x[3])*(x[2]+x[4]+2*G2*x[2]*x[4])+(-1+G2^2*x[1]*x[3])*(-1+G2^2*x[2]*x[4])
    F[2] = -x[3]-x[4]+x[2]*(-1+2*G1^2*x[3]*x[4]+G2^2*x[3]*x[4])+x[1]*(-1+G2^2*(x[3]*x[4]+x[2]*(x[3]+x[4]))+2*G1^2*(x[3]*x[4]+x[2]*(x[3]+x[4]+4*G2*x[3]*x[4])))
    F[3] = x[3]*x[4]+x[2]*(x[3]+x[4])+x[1]*(x[2]+x[3]+x[4]-4*G1^2*x[2]*x[3]*x[4]-2*G2^2*x[2]*x[3]*x[4])
    F[4] = 1/x[1]+1/x[2]+1/x[3]+1/x[4]
end
validsolution = "n"
while validsolution != "y"
    global s = nlsolve(f!, 10*(rand(nscatt)+rand(nscatt)*1im), ftol=1e-5, iterations=1_000)
    while s.f_converged == false
        s = nlsolve(f!, 10*(rand(nscatt)+rand(nscatt)*1im), ftol=1e-5, iterations=1_000)
    end
    s = nlsolve(f!, s.zero, ftol=1e-20, iterations=10_000)
    FF = ones(nscatt)+ones(nscatt)*1im
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

npzwrite("./data/eigs_4scatt_nlsolve.npz", Dict("rspan" => rspan, "rover2" => float(r), "epsilon" => ϵ, "alphas" => alphas, "scattpos" => scattpos, "eigs" => eigs))