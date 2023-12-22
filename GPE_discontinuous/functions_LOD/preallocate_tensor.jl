
function preallocate_tensor(V_LOD, M_LOD)

    N = size(V_LOD, 1)
    W_Iptr = Array{Int64}(undef, N + 1)
    W_J = Array{Int64}(undef, 5 * 10^7)
    W_K = Array{Int64}(undef, 5 * 10^7)

    W_Iptr[1] = 1
    N = size(M_LOD, 1)
    J = Array{Int64}(zeros(100 * N))
    K = similar(J)
    fill!(J, 0)
    Index = 0

    for basis = 1:N
        index = 0
        support = M_LOD[:, basis].nzind
        n_loc = length(support)
        for j = 1:n_loc
            if (support[j] >= basis)
                for k = j:n_loc
                    if (support[k] >= support[j])
                        if (abs(M_LOD[support[j], support[k]]) > 10^-16)
                            index += 1
                            J[index] = support[j]
                            K[index] = support[k]
                        end
                    end
                end
            end
        end

        W_J[Index+1:Index+index] = J[1:index]
        W_K[Index+1:Index+index] = K[1:index]
        Index += index
        W_Iptr[basis+1] = Index + 1

        if (mod(basis, 100) == 0)
            println(basis / N)
        end
        if (basis == N)
            println(Index)
        end
    end

    W_J = W_J[1:Index]
    W_K = W_K[1:Index]
    W_val = zeros(Index)

    return W_Iptr, W_J, W_K, W_val

end



