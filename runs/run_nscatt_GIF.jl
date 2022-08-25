include("../src/ScatterersDDA.jl")
using .ScatterersDDA
using LinearAlgebra
using NPZ
using ProgressMeter
using NLsolve

ω = 1.
J = Stdd(1.)
normalized = false
nscatt = 5
rsteps = 20

###################
### EP parameters
ωshift = 0.00001 # shift to lower the EP sensitivity and approach the inequality boundary to have a passive system
ϵ = 0.02
r = 2
rspan = LinRange(r-ϵ, r+ϵ, rsteps)

scattpos = zeros((rsteps, nscatt, 3))
for i in 1:rsteps
    scattpos[i, 1, :] = [rspan[i] 0. 0.]
    for j in 1:nscatt-1
        scattpos[i, j+1, :] = [r*cos(2π*j/nscatt) r*sin(2π*j/nscatt) 0.]
    end
end

# ϕinput = rand(nscatt).+rand(nscatt).*1im
ϕplane = r -> exp(1im*k0(ω, J)*r)
ϕinput  = zeros(ComplexF64, (rsteps, nscatt))
for k in 1:rsteps
    for i in 1:nscatt
        ϕinput[k, i] = ϕplane(scattpos[k, i, 1]) # plane wave propagating along the x-axis
    end
end

###################
### non-linear system solver
Gvect = zeros(floor(Int, nscatt/2))+zeros(floor(Int, nscatt/2))*1im
for i in eachindex(Gvect)
    Gvect[i] = greensfun(scattpos[1, 2, :], scattpos[1, 2+i, :], ω+ωshift*1im, J)
end
G1 = Gvect[1]
G2 = Gvect[2]
function f!(F, x)
    F[1] = 1 - 2*G1^5*x[1]*x[2]*x[3]*x[4]*x[5] - 2*G2^5*x[1]*x[2]*x[3]*x[4]*x[5] - G2^2*(x[1]*(x[3] + x[4]) + x[3]*x[5] + x[2]*(x[4] + x[5])) + G2^4*(x[2]*x[3]*x[4]*x[5] + x[1]*(x[2]*x[4]*x[5] + x[3]*x[4]*x[5] + x[2]*x[3]*(x[4] + x[5]))) - 2*G1^3*G2*(x[2]*x[3]*x[4]*x[5] + x[1]*(x[2]*x[4]*x[5] + x[3]*x[4]*x[5] + x[2]*x[3]*(x[4] + x[5] + 5*G2*x[4]*x[5]))) + G1^4*(x[2]*x[3]*x[4]*x[5] + x[1]*(x[2]*x[4]*x[5] + x[3]*x[4]*x[5] + x[2]*x[3]*(x[4] + x[5] + 10*G2*x[4]*x[5]))) + 2*G1*G2^2*((-x[2])*(x[3] + x[4] + G2*x[3]*x[4])*x[5] + x[1]*((-G2)*x[2]*x[3]*x[5] - x[3]*(x[4] + x[5] + G2*x[4]*x[5]) + x[2]*x[4]*(-1 + 5*G2^2*x[3]*x[5] - G2*(x[3] + x[5])))) - G1^2*(x[1]*(1 + 2*G2*x[4] + G2^2*x[3]*x[4])*x[5] + x[4]*(x[3] + x[5] + 2*G2*x[3]*x[5]) + x[2]*x[3]*(1 + 2*G2*x[4] + G2^2*x[4]*x[5]) + x[1]*x[2]*(1 + 10*G2^3*x[3]*x[4]*x[5] + 2*G2*(x[3] + x[5]) + G2^2*(x[4]*x[5] + x[3]*(x[4] + x[5]))))
    F[2] = -x[3] - x[4] - x[5] + 2*G1^2*x[3]*x[4]*x[5] + G2^2*x[3]*x[4]*x[5] + x[2]*(-1 + 4*G1*G2^2*x[3]*x[4]*x[5] + G2^2*(x[3]*x[4] + 2*x[3]*x[5] + 2*x[4]*x[5]) + G1^2*(x[4]*x[5] + x[3]*(2*x[4] + x[5] + 4*G2*x[4]*x[5]))) + x[1]*(-1 - 5*G1^4*x[2]*x[3]*x[4]*x[5] + 10*G1^3*G2*x[2]*x[3]*x[4]*x[5] - 5*G2^4*x[2]*x[3]*x[4]*x[5] + G2^2*(x[4]*x[5] + 2*x[3]*(x[4] + x[5]) + x[2]*(x[3] + 2*x[4] + x[5])) + 2*G1*G2^2*(2*x[3]*x[4]*x[5] + x[2]*(2*x[3]*x[5] + 2*x[4]*x[5] + x[3]*x[4]*(2 + 5*G2*x[5]))) + G1^2*(2*x[4]*x[5] + x[3]*(x[4] + x[5] + 4*G2*x[4]*x[5]) + x[2]*(x[4] + 2*x[5] + 4*G2*x[4]*x[5] + x[3]*(2 + 5*G2^2*x[4]*x[5] + 4*G2*(x[4] + x[5])))))
    F[3] = x[3]*x[4] + x[3]*x[5] + x[4]*x[5] + x[2]*(x[3] + x[4] + x[5] - 3*G1^2*x[3]*x[4]*x[5] - 3*G2^2*x[3]*x[4]*x[5]) + x[1]*(x[3] + x[4] + x[5] - 3*G1^2*x[3]*x[4]*x[5] - 3*G2^2*x[3]*x[4]*x[5] - x[2]*(-1 + 10*G1*G2^2*x[3]*x[4]*x[5] + 3*G2^2*(x[4]*x[5] + x[3]*(x[4] + x[5])) + G1^2*(3*x[3]*x[5] + 3*x[4]*x[5] + x[3]*x[4]*(3 + 10*G2*x[5]))))
    F[4] = x[1]*x[4]*x[5] + x[2]*x[4]*x[5] + x[3]*x[4]*x[5] + x[1]*x[3]*(x[4] + x[5]) + x[2]*x[3]*(x[4] + x[5]) + x[1]*x[2]*(x[3] + x[4] + x[5] - 5*G1^2*x[3]*x[4]*x[5] - 5*G2^2*x[3]*x[4]*x[5])
    F[5] = 1/x[1] + 1/x[2] + 1/x[3] + 1/x[4] + 1/x[5]
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
### Coalescence eigenvalues
# eigs = zeros(ComplexF64, (rsteps, length(alphas)))
# p = Progress(length(rspan));
# println("Starting...")
# Threads.@threads for k in eachindex(rspan)
#     eigs[k, :] = eigvals(intmatrix(scattpos[k, :, :], alphas, ω, J; normalized=normalized, imagshift=1E-23))
#     next!(p)
# end
# npzwrite("./data/eigs_n5_omegashift.npz", Dict("rspan" => rspan, "r" => float(r), "epsilon" => ϵ, "alphas" => alphas, "scattpos" => scattpos, "eigs" => eigs))

###################
### Total field over r
# xx = LinRange(-10, 10, 200)
# yy = LinRange(-10, 10, 200)
# phitot = zeros(ComplexF64, (rsteps, length(xx), length(yy)))
# p = Progress(rsteps);
# println("Starting...")
# Threads.@threads for k in eachindex(rspan)
#     for i in eachindex(xx)
#         for j in eachindex(yy)
#             phitot[k, j, i] = totfield([xx[i], yy[j], 0], ϕinput[k, :], scattpos[k, :, :], alphas, ω, J; normalized=normalized, imagshift=1E-23)
#         end
#     end
#     next!(p)
# end
# npzwrite("./data/totalfield_n5_omegashift.npz", Dict("rspan" => rspan, "r" => float(r), "epsilon" => ϵ, "alphas" => alphas, "scattpos" => scattpos, "phiinput" => ϕinput, "xx" => xx, "yy" => yy, "phitot" => phitot))

###################
### Power output over r
# maxradius = maximum([norm(scattpos[end, i, :]) for i in eachindex(scattpos[end, :, 1])])
maxradius = 10
ωsteps = 200
ωspan = LinRange(ω-0.01, ω+0.01, ωsteps)
Pout = zeros((rsteps, ωsteps))
p = Progress(rsteps);
println("Starting...")
Threads.@threads for k in 1:rsteps
    for i in 1:ωsteps
        Pout[k, i] = poweroutexplicit(maxradius, ϕinput[k, :], scattpos[k, :, :], alphas, ωspan[i], J; normalized=normalized)[1]
    end
    next!(p)
end
npzwrite("./data/poweroveromega_n5_omegashift.npz", Dict("rspan" => rspan, "omegaspan" => ωspan, "r" => float(r), "epsilon" => ϵ, "alphas" => alphas, "scattpos" => scattpos, "phiinput" => ϕinput, "Pout" => Pout))