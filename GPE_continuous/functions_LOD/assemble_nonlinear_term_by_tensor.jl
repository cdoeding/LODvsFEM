function assemble_density(W_Iptr, W_J, W_K, W_val, U, Rho)

    Uconj = U'

    for i = 1:length(W_Iptr)-1
        for idx = W_Iptr[i]:W_Iptr[i+1]-1
            j = W_J[idx]
            k = W_K[idx]

            if (i == k)
                Rho[i] += W_val[idx] * real(U[j] * Uconj[k])
            elseif (j == k)
                Rho[i] += W_val[idx] * real(U[j] * Uconj[k])
                Rho[j] += 2 * W_val[idx] * real(U[i] * Uconj[k])
            elseif (i == j)
                Rho[i] += 2 * W_val[idx] * real(U[j] * Uconj[k])
                Rho[k] += W_val[idx] * real(U[i] * Uconj[j])
            else
                Rho[i] += 2 * W_val[idx] * real(U[j] * Uconj[k])
                Rho[j] += 2 * W_val[idx] * real(U[i] * Uconj[k])
                Rho[k] += 2 * W_val[idx] * real(U[i] * Uconj[j])
            end
        end
    end

    return nothing

end

function assemble_nonlinear_term(W_Iptr, W_J, W_K, W_val, U, Rho, NL)

    for i = 1:length(W_Iptr)-1
        for idx = W_Iptr[i]:W_Iptr[i+1]-1
            j = W_J[idx]
            k = W_K[idx]

            if (i == k)
                NL[i] += W_val[idx] * (Rho[j] * U[k])
            elseif (j == k)
                NL[i] += W_val[idx] * (Rho[j] * U[k])
                NL[j] += W_val[idx] * (Rho[i] * U[k] + Rho[k] * U[i])
            elseif (i == j)
                NL[i] += W_val[idx] * (Rho[j] * U[k] + Rho[k] * U[j])
                NL[k] += W_val[idx] * (Rho[i] * U[j])
            else
                NL[i] += W_val[idx] * (Rho[j] * U[k] + Rho[k] * U[j])
                NL[j] += W_val[idx] * (Rho[i] * U[k] + Rho[k] * U[i])
                NL[k] += W_val[idx] * (Rho[i] * U[j] + Rho[j] * U[i])
            end
        end
    end

    return nothing
end



