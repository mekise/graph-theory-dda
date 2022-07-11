function integrand(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    (r, θ, ϕ) = (x[1], x[2], x[3])
    grad = complexgrad(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)
    arg = conj(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))*[1, 1/r, 1/r/sin(θ)].*grad
    return imag(arg)*r^2*sin(θ)
end

function powerout(r, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    f(ϕ, θ) = integrand([r, θ, ϕ], ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)[1]
    return dblquad(f, 0, 2π, θ->0, θ->π)
end