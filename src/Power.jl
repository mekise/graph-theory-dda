function integrand(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    grad = x -> complexder(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)
    arg = x -> conj(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))*grad(x)
    return imag(arg(x))*(x[1]^2)*sin(x[2])
end

# function powerout(r, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
#     f(ϕ, θ) = integrand([r, θ, ϕ], ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)[1]
#     return dblquad(f, 0, 2π, θ->0, θ->π)
# end

function powerout(r, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    f(x) = integrand([r, x[1], x[2]], ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)
    return hcubature(f, [0, 0], [π, 2π])
end

#########################

function greensfunpolar(a, ω, J::Stdd; imagshift=1E-12)
    (r, θ, ϕ) = (a[1], a[2], a[3])
    return greensfun([r*cos(θ)*sin(ϕ), r*sin(θ)*sin(ϕ), r*cos(ϕ)], [0., 0., 0.], ω, J; imagshift=imagshift)
end

function Gcomplexgrad(x, ω, J::Stdd; imagshift=1E-12)
    freal = x -> real(greensfunpolar(x, ω, J; imagshift=imagshift))
    fimag = x -> imag(greensfunpolar(x, ω, J; imagshift=imagshift))
    gradreal = x -> ForwardDiff.gradient(freal, x)
    gradimag = x -> ForwardDiff.gradient(fimag, x)
    return gradreal(x).+1im*gradimag(x)
end

function Gintegrand(x, ω, J::Stdd; imagshift=1E-12)
    (r, θ) = (x[1], x[2])
    return imag(conj(greensfunpolar(x, ω, J; imagshift=imagshift))*Gcomplexgrad(x, ω, J; imagshift=imagshift))*(r^2)*sin(θ)
end

function Gintegral(r, ω, J::Stdd; imagshift=1E-12)
    f(ϕ, θ) = Gintegrand([r, θ, ϕ], ω, J; imagshift=imagshift)[1]
    return dblquad(f, 0, 2π, θ->0, θ->π)
end