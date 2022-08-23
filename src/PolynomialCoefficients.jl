function buildCoefficients(scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-23)
    n = size(scattpos, 1)
    M = intmatrix(scattpos, alphas, ω, J; normalized=normalized, imagshift=imagshift)

    res = λ
end