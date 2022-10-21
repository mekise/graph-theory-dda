function intmatrix(scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    n = size(scattpos, 1)
    M = zeros(ComplexF64, (n, n))
    for i in 1:n
        for j in 1:n
            M[i, j] = -greensfun(scattpos[i, :], scattpos[j, :], ω, J; imagshift=imagshift)
        end
    end
    M[diagind(M)] .= 1 ./ alphas
    return M
end

function weakapproxmatrix(scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    M = intmatrix(scattpos, x, ω, J)
    det(M[setdiff(1:end, i), setdiff(1:end, j)])
    sol = sum(det.([intmatrix(scattpos, x, ω, J)[setdiff(1:end, i), setdiff(1:end, j)] for i in combinations([j for j in 1:nscatt], cindex-1)]))


    M-1 = C/det


end