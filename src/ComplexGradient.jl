function complexgrad(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-23)
    freal = x -> real(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))
    fimag = x -> imag(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))
    gradreal = x -> ForwardDiff.gradient(freal, x)
    gradimag = x -> ForwardDiff.gradient(fimag, x)
    return gradreal(x).+1im*gradimag(x)
end

function complexder(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-23)
    r, θ, ϕ = x[1], x[2], x[3]
    freal = r -> real(totfieldpolar([r, θ, ϕ], ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))
    fimag = r -> imag(totfieldpolar([r, θ, ϕ], ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))
    gradreal = r -> ForwardDiff.derivative(freal, r)
    gradimag = r -> ForwardDiff.derivative(fimag, r)
    return gradreal(r).+1im*gradimag(r)
end