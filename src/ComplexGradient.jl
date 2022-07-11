function complexgrad(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    freal = x -> real(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))
    fimag = x -> imag(totfieldpolar(x, ϕinput, scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift))
    gradreal = x -> ForwardDiff.gradient(freal, x)
    gradimag = x -> ForwardDiff.gradient(fimag, x)
    return gradreal(x).+1im*gradimag(x)
end