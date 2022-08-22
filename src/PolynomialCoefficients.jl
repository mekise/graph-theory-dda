function buildCoefficients(scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-23)
    n = size(scattpos, 1)
    M = intmatrix(scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)

    res = λ
end

G1 = greensfun(scattpos[1, 1, :], scattpos[1, 5, :], ω, J)
G2 = greensfun(scattpos[1, 1, :], scattpos[1, 4, :], ω, J)
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
