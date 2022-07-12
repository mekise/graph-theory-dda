function integrand(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    grad = x -> complexder(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)
    arg = x -> conj(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))*grad(x)
    return imag(arg(x))*(x[1]^2)*sin(x[2])
end

function powerout(r, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    f(x) = integrand([r, x[1], x[2]], ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)
    return hcubature(f, [0, 0], [π, 2π]; reltol=1E-7)
end

function evalsumm(nscatt, ϕinput, scattpos, alphas, ω, J::Stdd)
    summ = 0
    for i in 1:nscatt
        summ += abs(incfield(ϕinput, scattpos, alphas, ω, J)[i])^2 * imag(alphas[i])
    end
    return -summ
end