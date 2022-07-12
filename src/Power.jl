function integrand(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-23)
    grad = x -> complexder(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)
    arg = x -> conj(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))*grad(x)
    return imag(arg(x))
end

function powerout(r, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-23, reltol=1e-8)
    f(x) = integrand([r, x[1], x[2]], ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)*(r^2)*sin(x[1])
    return hcubature(f, [0, 0], [π, 2π]; reltol=reltol)
end

function evalsumm(nscatt, ϕinput, scattpos, alphas, ω, J::Stdd)
    summ = 0
    ϕinc = incfield(ϕinput, scattpos, alphas, ω, J)
    for i in 1:nscatt
        summ += abs(ϕinc[i])^2 * imag(alphas[i])
    end
    return -summ
end