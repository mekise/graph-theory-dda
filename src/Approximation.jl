function determinant(M)
    indices = [i for i in 1:size(M)[1]]
    result = 0. + 0im
    if size(M) == (2,2)
        return (M[1,1]*M[2,2]-M[2,1]*M[1,2])
    end
    for i in indices
        result += (-1)^(i+1) * M[1, i] * determinant(M[setdiff(1:end,1), setdiff(1:end,i)])
    end
    return result
end

function det_weak(M)
    n = size(M)[1]
    indices = [i for i in 1:n]
    result = 0. + 0im
    result += prod(M[diagind(M)])
    perms = permsingleinversion(M)
    perms = [perms[:,i] for i in axes(perms, 2)]
    for i in perms
        product = 1. + 0im 
        for k in 1:n
            product *= M[i[k], indices[k]]
        end
        result += -product
    end
    return result
end

function permsingleinversion(n::Int)
    indices = [i for i in 1:n]
    result = zeros(Int64, 0)
    count = 0
    for i in 1:n
        for j in i+1:n
            tmp = copy(indices)
            tmp[i], tmp[j] = tmp[j], tmp[i]
            append!(result, tmp)
            count += 1
        end
    end
    return reshape(result, (n, count))
    # return [result[:, i] for i in axes(result, 2)]
end

function permsingleinversion(M)
    n = size(M)[1]
    indices = [i for i in 1:n]
    result = zeros(Int64, 0)
    count = 0
    for i in 1:n
        for j in i+1:n
            tmp = copy(indices)
            tmp[i], tmp[j] = tmp[j], tmp[i]
            append!(result, tmp)
            count += 1
        end
    end
    return reshape(result, (n, count))
    # return [result[:, i] for i in axes(result, 2)]
end

function weak4by4(scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    G1 = greensfun(scattpos[1,:], scattpos[2,:], ω, J; imagshift=imagshift)
    G2 = greensfun(scattpos[1,:], scattpos[3,:], ω, J; imagshift=imagshift)
    matrix = zeros(ComplexF64, (4,4))
    matrix[1,1] = -(G1^2/alphas[1]) - G2^2/alphas[2] - G1^2/alphas[3] + 1/(alphas[1]*alphas[2]*alphas[3])
    matrix[1,2] = -1 * (-((G1*G2)/alphas[1])-(G1*G2)/alphas[2]-G1/(alphas[1]*alphas[2]))
    matrix[1,3] = G1^2/alphas[1]+G1^2/alphas[3]+G2/(alphas[1]*alphas[3])
    matrix[1,4] = -1 * (-((G1*G2)/alphas[2]) - (G1*G2)/alphas[3] - G1/(alphas[2]*alphas[3]))
    matrix[2,2] = -(G2^2/alphas[1]) - G1^2/alphas[2] - G1^2/alphas[4] + 1/(alphas[1]*alphas[2]*alphas[4])
    matrix[2,3] = -1 * (-((G1*G2)/alphas[1]) - (G1*G2)/alphas[4] - G1/(alphas[1]*alphas[4]))
    matrix[2,4] = +(G1^2/alphas[2]) + G1^2/alphas[4] + G2/(alphas[2]*alphas[4])
    matrix[3,3] = -(G1^2/alphas[1]) - G1^2/alphas[3] - G2^2/alphas[4] + 1/(alphas[1]*alphas[3]*alphas[4])
    matrix[3,4] = -1 * (-((G1*G2)/alphas[3]) - (G1*G2)/alphas[4] - G1/(alphas[3]*alphas[4]))
    matrix[4,4] = -(G1^2/alphas[2]) - G2^2/alphas[3] - G1^2/alphas[4] + 1/(alphas[2]*alphas[3]*alphas[4])
    for i in 1:4
        for j in i+1:4
            matrix[j,i] = matrix[i,j]
        end
    end
    det = -(G1^2/(alphas[1]*alphas[2])) - G2^2/(alphas[1]*alphas[3]) - G1^2/(alphas[2]*alphas[3]) - G1^2/(alphas[1]*alphas[4]) - G2^2/(alphas[2]*alphas[4]) - G1^2/(alphas[3]*alphas[4]) + 1/(alphas[1]*alphas[2]*alphas[3]*alphas[4])
    matrix ./= det
end

function weakincfield(ϕinput, scattpos, alphas, ω, J::Stdd; returnfield=false, imagshift=1E-23)
    M = weak4by4(scattpos, alphas, ω, J; imagshift=imagshift)
    ϕinc_tilde = similar(ϕinput)
    mul!(ϕinc_tilde, M, (ϕinput./alphas))
    if returnfield
        return ϕinc_tilde.*alphas
    else
        return ϕinc_tilde
    end
end

function weaktotfield(x, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    ϕinc_tilde = weakincfield(ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)
    n = size(scattpos, 1)
    ϕtot = 0. + 0im
    for i in 1:n
        ϕtot += ϕinc_tilde[i]*greensfun(x, scattpos[i, :], ω, J; imagshift=imagshift)
    end
    return ϕtot
end

function strong4by4(scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    G1 = greensfun(scattpos[1,:], scattpos[2,:], ω, J; imagshift=imagshift)
    G2 = greensfun(scattpos[1,:], scattpos[3,:], ω, J; imagshift=imagshift)
    matrix = zeros(ComplexF64, (4,4))
    matrix[1,1] = -2*G1^2*G2
    matrix[1,2] = G1*G2^2
    matrix[1,3] = 2*G1^2*G2-G2^3
    matrix[1,4] = G1*G2^2
    matrix[2,2] = -2*G1^2*G2
    matrix[2,3] = G1*G2^2
    matrix[2,4] = 2*G1^2*G2-G2^3
    matrix[3,3] = -2*G1^2*G2
    matrix[3,4] = G1*G2^2
    matrix[4,4] = -2*G1^2*G2
    for i in 1:4
        for j in i+1:4
            matrix[j,i] = matrix[i,j]
        end
    end
    det = G2^4-4*G1^2*G2^2-2*G1^2*G2*(1/alphas[1]+1/alphas[2]+1/alphas[3]+1/alphas[4])
    matrix ./= det
end

function strongincfield(ϕinput, scattpos, alphas, ω, J::Stdd; returnfield=false, imagshift=1E-23)
    M = strong4by4(scattpos, alphas, ω, J; imagshift=imagshift)
    ϕinc_tilde = similar(ϕinput)
    mul!(ϕinc_tilde, M, (ϕinput./alphas))
    if returnfield
        return ϕinc_tilde.*alphas
    else
        return ϕinc_tilde
    end
end

function strongtotfield(x, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    ϕinc_tilde = strongincfield(ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)
    n = size(scattpos, 1)
    ϕtot = 0. + 0im
    for i in 1:n
        ϕtot += ϕinc_tilde[i]*greensfun(x, scattpos[i, :], ω, J; imagshift=imagshift)
    end
    return ϕtot
end

# function cofactor(a, b, scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
#     G1 = greensfun(scattpos[1], scattpos[2], ω, J; imagshift=imagshift)
#     G2 = greensfun(scattpos[1], scattpos[3], ω, J; imagshift=imagshift)
#     det = prod(alphas.^(-1)) - (alphas[1]^(-1)*alphas[2]^(-1)+alphas[2]^(-1)*alphas[3]^(-1)+alphas[3]^(-1)*alphas[4]^(-1)+alphas[4]^(-1)*alphas[1]^(-1))*G1^2 - (alphas[1]^(-1)*alphas[3]^(-1)+alphas[2]^(-1)*alphas[4]^(-1))*G2^2
#     if a==b
#         tmp_alphas = alphas[setdiff(1:end, a)]
#         if a==1
#             return (prod(tmp_alphas.^(-1)) - (alphas[end]^(-1)+alphas[a+1])*G1^2 - (alphas[a+2]^2)*G2^2)
#         elseif a==length(alphas)
#             return (prod(tmp_alphas.^(-1)) - (alphas[a-1]^(-1)+alphas[1])*G1^2 - (alphas[a+2]^2)*G2^2)
#         else
#             return (prod(tmp_alphas.^(-1)) - (alphas[a-1]^(-1)+alphas[a+1])*G1^2 - (alphas[a+2]^2)*G2^2)
#         end
#     else
#         return  (-1)^(a+b) * ()
#     end
# end

# function intmatrix_weak(scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
#     nscatt = length(alphas)
#     result = zeros(ComplexF64, (nscatt, nscatt))
#     for i in 1:nscatt
#         for j in 1:nscatt
#             result[i, j] = cofactor(i, j, scattpos, alphas, ω, J; imagshift=imagshift)
#         end
#     end
#     return result
# end