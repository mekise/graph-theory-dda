function integrand(x, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    grad = x -> complexder(x, ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)
    arg = x -> conj(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; imagshift=imagshift))*grad(x)
    return imag(arg(x))
end

function powerout(r, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23, reltol=1e-8)
    f(x) = integrand([r, x[1], x[2]], ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)*(r^2)*sin(x[1])
    return hcubature(f, [0, 0], [π, 2π]; reltol=reltol)
end

function poweroutexplicit(r, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23, reltol=1e-8)
    summ = 0
    incf = incfield(ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)
    for i in 1:length(scattpos[:, 1])
        for j in 1:length(scattpos[:, 1])
            function f(x)
                ri = norm([r*cos(x[1])*sin(x[2]), r*sin(x[1])*sin(x[2]), r*cos(x[2])] - scattpos[i, :])
                rj = norm([r*cos(x[1])*sin(x[2]), r*sin(x[1])*sin(x[2]), r*cos(x[2])] - scattpos[j, :])
                prefactor = conj(alphas[i])*alphas[j]*conj(incf[i])*incf[j]
                return imag(prefactor * (exp(1im*(rj-ri+π/2))/(rj*ri) - exp(1im*(ri-rj))/(ri*rj^2))) *(r^2)*sin(x[1])
            end
            integral = hcubature(f, [0, 0], [π, 2π]; reltol=reltol)[1]
            summ += integral /(4π)^2
        end
    end
    return summ
end

function evalsumm(nscatt, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    summ = 0
    ϕinc = incfield(ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)
    for i in 1:nscatt
        summ += - abs(ϕinc[i])^2 * imag(alphas[i])
    end
    return summ
end

function evalsummcorrected(nscatt, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    summ = 0
    ϕinc = incfield(ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)
    for i in 1:nscatt
        summ += - abs(ϕinc[i])^2 * imag(alphas[i]) + imag(conj(ϕinc[i]) * ϕinput[i])
    end
    return summ
end