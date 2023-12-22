function solve_GPE_1D_q2(u0, k, Nt, A_LOD, M_LOD, M_V_LOD, W_Iptr, W_J, W_K, W_val, alpha, beta, tol, max_it)

    # cG(q)-scheme parameter
    q = 2 #degree in time
    tau, w = get_gauss_nodes(2, 0, 1) #Gauss integration
    tau_t, w_t = get_gauss_nodes(2 * q, 0, 1) # Gauss integration for nonlinear term

    # Lagrange polynomials
    li = [x -> (x - tau[2]) ./ (tau[1] - tau[2]), x -> (x - tau[1]) ./ (tau[2] - tau[1])]
    li_hat = [x -> (x - tau[1]) .* (x - tau[2]) ./ (tau[1] * tau[2]);
        x -> x .* (x - tau[2]) ./ (tau[1]^2 - tau[1] * tau[2]);
        x -> x .* (x - tau[1]) ./ (tau[2]^2 - tau[1] * tau[2])]

    # Derivatives of Lagrange polynomials
    dli = [x -> 1 / (tau[1] - tau[2]); x -> 1 / (tau[2] - tau[1])]

    # compute RK coefficients
    m = zeros(q, q)
    for i = 1:q
        for j = 1:q
            if i == j
                m[i, j] = (w[i] / tau[j]) * (1 + tau[i] * dli[j](tau[i]))
            else
                m[i, j] = (w[i] / tau[j]) * tau[i] * dli[j](tau[i])
            end
        end
    end

    mi = m \ I
    gamma, Si = eigen(mi * diagm(w))
    S = Si \ I

    # a coefficients
    a = zeros(q, 1) * 0im
    for i = 1:q
        a[i] = sum(S[i, :])
    end

    # b coeeficitents
    b = zeros(q, 2 * q) * 0im
    SMinv = S * mi
    for i = 1:q
        for m = 1:2*q
            for j = 1:q
                b[i, m] = b[i, m] + SMinv[i, j] * li[j](tau_t[m])
            end
            b[i, m] = w_t[m] * b[i, m]
        end
    end
    # c coefficients
    c = zeros(q, 2 * q) * 0im
    for j = 1:q
        for m = 1:2*q
            for i = 1:q
                c[j, m] = c[j, m] + li_hat[i+1](tau_t[m]) * Si[i, j]
            end
        end
    end
    c0 = zeros(2 * q, 1)
    for m = 1:2*q
        c0[m] = li_hat[1](tau_t[m])
    end

    # compute LHS and LU-decomposition
    LU_Matrix_1 = lu(M_LOD + alpha * k * gamma[1] * A_LOD + 1im * k * gamma[1] * M_V_LOD)
    LU_Matrix_2 = lu(M_LOD + alpha * k * gamma[2] * A_LOD + 1im * k * gamma[2] * M_V_LOD)

    # time integration
    # set intial starting point for fixed point iteration
    dim = size(M_LOD, 1)
    xi_old_1 = zeros(dim) * 0im
    xi_old_2 = zeros(dim) * 0im


    rho_interim = zeros(dim)
    rho_LOD_interim = zeros(dim)

    f1 = zeros(dim) * 0im
    f2 = zeros(dim) * 0im
    f3 = zeros(dim) * 0im
    f4 = zeros(dim) * 0im

    U_sol = zeros(dim, Nt + 1) * 0im
    U = zeros(dim) * 0im
    copyto!(U, u0)
    U_sol[:, 1] = U

    Mlu = lu(M_LOD)

    for n = 1:Nt

        # fixed-point itreration
        delta = 1
        l = 0

        copyto!(xi_old_1, U)
        copyto!(xi_old_2, U)

        while delta > tol && l < max_it

            U_interim = c0[1] * U + c[1, 1] * xi_old_1 + c[2, 1] * xi_old_2
            fill!(rho_interim, 0.0)
            assemble_density(W_Iptr, W_J, W_K, W_val, U_interim, rho_interim)
            rho_LOD_interim = Mlu \ rho_interim
            fill!(f1, 0.0)
            assemble_nonlinear_term(W_Iptr, W_J, W_K, W_val, U_interim, rho_LOD_interim, f1)

            U_interim = c0[2] * U + c[1, 2] * xi_old_1 + c[2, 2] * xi_old_2
            fill!(rho_interim, 0.0)
            assemble_density(W_Iptr, W_J, W_K, W_val, U_interim, rho_interim)
            rho_LOD_interim = Mlu \ rho_interim
            fill!(f2, 0.0)
            assemble_nonlinear_term(W_Iptr, W_J, W_K, W_val, U_interim, rho_LOD_interim, f2)

            U_interim = c0[3] * U + c[1, 3] * xi_old_1 + c[2, 3] * xi_old_2
            fill!(rho_interim, 0.0)
            assemble_density(W_Iptr, W_J, W_K, W_val, U_interim, rho_interim)
            rho_LOD_interim = Mlu \ rho_interim
            fill!(f3, 0.0)
            assemble_nonlinear_term(W_Iptr, W_J, W_K, W_val, U_interim, rho_LOD_interim, f3)

            U_interim = c0[4] * U + c[1, 4] * xi_old_1 + c[2, 4] * xi_old_2
            fill!(rho_interim, 0.0)
            assemble_density(W_Iptr, W_J, W_K, W_val, U_interim, rho_interim)
            rho_LOD_interim = Mlu \ rho_interim
            fill!(f4, 0.0)
            assemble_nonlinear_term(W_Iptr, W_J, W_K, W_val, U_interim, rho_LOD_interim, f4)

            # compute right hand side
            F1 = a[1] * M_LOD * U - 1im * k * beta * (b[1, 1] * f1 + b[1, 2] * f2 + b[1, 3] * f3 + b[1, 4] * f4)
            F2 = a[2] * M_LOD * U - 1im * k * beta * (b[2, 1] * f1 + b[2, 2] * f2 + b[2, 3] * f3 + b[2, 4] * f4)

            xi_new_1 = LU_Matrix_1 \ F1
            xi_new_2 = LU_Matrix_2 \ F2

            # update iteration paramter
            dU_1 = xi_old_1 - xi_new_1
            dU_2 = xi_old_2 - xi_new_2
            delta = sqrt(max(real(dot(dU_1, M_LOD * dU_1)), real(dot(dU_2, M_LOD * dU_2))))

            copyto!(xi_old_1, xi_new_1)
            copyto!(xi_old_2, xi_new_2)

            l += 1

            print("->")
            print(l)

        end

        xi_transform_1 = Si[1, 1] * xi_old_1 + Si[1, 2] * xi_old_2
        xi_transform_2 = Si[2, 1] * xi_old_1 + Si[2, 2] * xi_old_2
        U = li_hat[1](1) * U + li_hat[2](1) * xi_transform_1 + li_hat[3](1) * xi_transform_2
        U_sol[:, n+1] = U
        println(" ")
        println(n)
    end

    return U_sol

end

