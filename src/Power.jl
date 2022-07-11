function integrand(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    (r, θ, ϕ) = (x[1], x[2], x[3])
    grad = complexgrad(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)
    arg = conj(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))*[1, 1/r, 1/r/sin(θ)].*grad
    return imag(arg)*(r^2)*sin(θ)
end

function powerout(r, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    f(ϕ, θ) = integrand([r, θ, ϕ], ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)[1]
    return dblquad(f, 0, 2π, θ->0, θ->π)
end

#########################

function greensfunpolar(a, ω, J::Stdd; imagshift=1E-12)
    (r, θ, ϕ) = (a[1], a[2], a[3])
    return greensfun([r*cos(θ)*sin(ϕ), r*sin(θ)*sin(ϕ), r*cos(ϕ)], [0., 0., 0.], ω, J; imagshift=imagshift)
end

function Gcomplexgrad(r, ω, J::Stdd; imagshift=1E-12)
    freal = r -> real(greensfunpolar([r, 0, 0], ω, J; imagshift=imagshift))
    fimag = r -> imag(greensfunpolar([r, 0, 0], ω, J; imagshift=imagshift))
    gradreal = r -> ForwardDiff.gradient(freal, r)
    gradimag = r -> ForwardDiff.gradient(fimag, r)
    return gradreal(r).+1im*gradimag(r)
end

function Gintegrand(r, ω, J::Stdd; imagshift=1E-12)
    return imag(conj(greensfunpolar(r, ω, J; imagshift=imagshift))*Gcomplexgrad(r, ω, J; imagshift=imagshift))*(x[1]^2)*sin(x[2])
end

function Gintegral(r, ω, J::Stdd; imagshift=1E-12)
    f(ϕ, θ) = Gintegrand([r, θ, ϕ], ω, J; imagshift=imagshift)[1]
    return dblquad(f, 0, 2π, θ->0, θ->π)
end