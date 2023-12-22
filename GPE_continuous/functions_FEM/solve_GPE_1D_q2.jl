function solve_GPE_1D_q2(A, M, M_V, alpha, beta, V, u0, T, Nt, Tri, nodes, nodes2mesh, threshold, l_max)

    # time discretization
    k = T / Nt
    t = Array(LinRange(0, T, Nt + 1))

    # cG(q)-scheme parameter
    q = 2 #degree in time
    tau, w = get_gauss_nodes(q, 0, 1) #Gauss integration
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
    D = diagm(gamma)
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

    # memory preallocation
    U = zeros(Nx - 1, Nt + 1) * im

    #-------------------- calculation
    # set intial value
    println(size(U))
    println(size(u0))
    U[:, 1] .= u0

    # compute LHS and LU-decomposition
    LU_Matrix_1 = lu(M + alpha * k * gamma[1] * A + 1im * k * gamma[1] * M_V)
    LU_Matrix_2 = lu(M + alpha * k * gamma[2] * A + 1im * k * gamma[2] * M_V)

    # time integration
    xi_old_1 = zeros(Nx - 1) * 0im
    xi_old_2 = zeros(Nx - 1) * 0im

    f1 = zeros(Nx - 1) * 0im
    f2 = zeros(Nx - 1) * 0im
    f3 = zeros(Nx - 1) * 0im
    f4 = zeros(Nx - 1) * 0im

    # set intial starting point for fixed point iteration
    copyto!(xi_old_1, U[:, 1])
    copyto!(xi_old_2, U[:, 1])

    for n = 1:Nt

        # fixed-point itreration
        delta = 1
        l = 0
        while delta > threshold && l < l_max

            # assemble F(xi_old)
            f1 = assemble_F(c0[1] * U[:, n] + c[1, 1] * xi_old_1 + c[2, 1] * xi_old_2, Tri, nodes, nodes2mesh)
            f2 = assemble_F(c0[2] * U[:, n] + c[1, 2] * xi_old_1 + c[2, 2] * xi_old_2, Tri, nodes, nodes2mesh)
            f3 = assemble_F(c0[3] * U[:, n] + c[1, 3] * xi_old_1 + c[2, 3] * xi_old_2, Tri, nodes, nodes2mesh)
            f4 = assemble_F(c0[4] * U[:, n] + c[1, 4] * xi_old_1 + c[2, 4] * xi_old_2, Tri, nodes, nodes2mesh)

            # compute right hand side
            F1 = a[1] * M * U[:, n] - 1im * k * beta * (b[1, 1] * f1 + b[1, 2] * f2 + b[1, 3] * f3 + b[1, 4] * f4)
            F2 = a[2] * M * U[:, n] - 1im * k * beta * (b[2, 1] * f1 + b[2, 2] * f2 + b[2, 3] * f3 + b[2, 4] * f4)


            # solve linear equations
            xi_new_1 = LU_Matrix_1 \ F1
            xi_new_2 = LU_Matrix_2 \ F2

            # update iteration paramter
            dU_1 = xi_old_1 - xi_new_1
            dU_2 = xi_old_2 - xi_new_2
            delta = sqrt(real(dot(dU_1, M * dU_1) + dot(dU_2, M * dU_2)))
            copyto!(xi_old_1, xi_new_1)
            copyto!(xi_old_2, xi_new_2)

            l = l + 1

            print("->")
            print(l)

        end

        # error message and break if iteration didnt converge
        if l == l_max
            println("no convergence")
            break
        end

        if isnan(delta)
            println("no convergence")
            break
        end

        # set result and calculate U(n+1)
        xi_transform_1 = Si[1, 1] * xi_old_1 + Si[1, 2] * xi_old_2
        xi_transform_2 = Si[2, 1] * xi_old_1 + Si[2, 2] * xi_old_2

        U[:, n+1] .= li_hat[1](1) * U[:, n] + li_hat[2](1) * xi_transform_1 + li_hat[3](1) * xi_transform_2

        # intermediate output
        println(" ")
        println(n)
    end

    return U

end

function solve_GPE_1D_q2_P2(A, M, M_V, alpha, beta, V, u0, T, Nt, Tri, nodes, nodes2mesh, threshold, l_max)

    # time discretization
    k = T / Nt
    t = Array(LinRange(0, T, Nt + 1))

    # cG(q)-scheme parameter
    q = 2 #degree in time
    tau, w = get_gauss_nodes(q, 0, 1) #Gauss integration
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
    D = diagm(gamma)
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

    # memory preallocation
    dim = size(A, 1)
    U = zeros(dim, Nt + 1) * im

    #-------------------- calculation
    # set intial value
    U[:, 1] .= u0

    # compute LHS and LU-decomposition
    LU_Matrix_1 = lu(M + alpha * k * gamma[1] * A + 1im * k * gamma[1] * M_V)
    LU_Matrix_2 = lu(M + alpha * k * gamma[2] * A + 1im * k * gamma[2] * M_V)

    # time integration
    xi_old_1 = zeros(dim) * 0im
    xi_old_2 = zeros(dim) * 0im

    f1 = zeros(dim) * 0im
    f2 = zeros(dim) * 0im
    f3 = zeros(dim) * 0im
    f4 = zeros(dim) * 0im

    # set intial starting point for fixed point iteration
    copyto!(xi_old_1, U[:, 1])
    copyto!(xi_old_2, U[:, 1])

    for n = 1:Nt

        # fixed-point itreration
        delta = 1
        l = 0
        while delta > threshold && l < l_max

            # assemble F(xi_old)
            f1 = assemble_F_P2(c0[1] * U[:, n] + c[1, 1] * xi_old_1 + c[2, 1] * xi_old_2, Tri, nodes, nodes2mesh)
            f2 = assemble_F_P2(c0[2] * U[:, n] + c[1, 2] * xi_old_1 + c[2, 2] * xi_old_2, Tri, nodes, nodes2mesh)
            f3 = assemble_F_P2(c0[3] * U[:, n] + c[1, 3] * xi_old_1 + c[2, 3] * xi_old_2, Tri, nodes, nodes2mesh)
            f4 = assemble_F_P2(c0[4] * U[:, n] + c[1, 4] * xi_old_1 + c[2, 4] * xi_old_2, Tri, nodes, nodes2mesh)

            # compute right hand side
            F1 = a[1] * M * U[:, n] - 1im * k * beta * (b[1, 1] * f1 + b[1, 2] * f2 + b[1, 3] * f3 + b[1, 4] * f4)
            F2 = a[2] * M * U[:, n] - 1im * k * beta * (b[2, 1] * f1 + b[2, 2] * f2 + b[2, 3] * f3 + b[2, 4] * f4)


            # solve linear equations
            xi_new_1 = LU_Matrix_1 \ F1
            xi_new_2 = LU_Matrix_2 \ F2

            # update iteration paramter
            dU_1 = xi_old_1 - xi_new_1
            dU_2 = xi_old_2 - xi_new_2
            delta = sqrt(real(dot(dU_1, M * dU_1) + dot(dU_2, M * dU_2)))
            copyto!(xi_old_1, xi_new_1)
            copyto!(xi_old_2, xi_new_2)

            l = l + 1

            print("->")
            print(l)

        end

        # error message and break if iteration didnt converge
        if l == l_max
            println("no convergence")
            break
        end

        if isnan(delta)
            println("no convergence")
            break
        end

        # set result and calculate U(n+1)
        xi_transform_1 = Si[1, 1] * xi_old_1 + Si[1, 2] * xi_old_2
        xi_transform_2 = Si[2, 1] * xi_old_1 + Si[2, 2] * xi_old_2

        U[:, n+1] .= li_hat[1](1) * U[:, n] + li_hat[2](1) * xi_transform_1 + li_hat[3](1) * xi_transform_2

        # intermediate output
        println(" ")
        println(n)
    end

    return U

end

function solve_GPE_1D_q2_P3(A, M, M_V, alpha, beta, V, u0, T, Nt, Tri, nodes, nodes2mesh, threshold, l_max)

    # time discretization
    k = T / Nt
    t = Array(LinRange(0, T, Nt + 1))

    # cG(q)-scheme parameter
    q = 2 #degree in time
    tau, w = get_gauss_nodes(q, 0, 1) #Gauss integration
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
    D = diagm(gamma)
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

    # memory preallocation
    dim = size(A, 1)
    U = zeros(dim, Nt + 1) * im

    #-------------------- calculation
    # set intial value
    U[:, 1] .= u0

    # compute LHS and LU-decomposition
    LU_Matrix_1 = lu(M + alpha * k * gamma[1] * A + 1im * k * gamma[1] * M_V)
    LU_Matrix_2 = lu(M + alpha * k * gamma[2] * A + 1im * k * gamma[2] * M_V)

    # time integration
    xi_old_1 = zeros(dim) * 0im
    xi_old_2 = zeros(dim) * 0im

    f1 = zeros(dim) * 0im
    f2 = zeros(dim) * 0im
    f3 = zeros(dim) * 0im
    f4 = zeros(dim) * 0im

    # set intial starting point for fixed point iteration
    copyto!(xi_old_1, U[:, 1])
    copyto!(xi_old_2, U[:, 1])

    for n = 1:Nt

        # fixed-point itreration
        delta = 1
        l = 0
        while delta > threshold && l < l_max

            # assemble F(xi_old)
            f1 = assemble_F_P3(c0[1] * U[:, n] + c[1, 1] * xi_old_1 + c[2, 1] * xi_old_2, Tri, nodes, nodes2mesh)
            f2 = assemble_F_P3(c0[2] * U[:, n] + c[1, 2] * xi_old_1 + c[2, 2] * xi_old_2, Tri, nodes, nodes2mesh)
            f3 = assemble_F_P3(c0[3] * U[:, n] + c[1, 3] * xi_old_1 + c[2, 3] * xi_old_2, Tri, nodes, nodes2mesh)
            f4 = assemble_F_P3(c0[4] * U[:, n] + c[1, 4] * xi_old_1 + c[2, 4] * xi_old_2, Tri, nodes, nodes2mesh)

            # compute right hand side
            F1 = a[1] * M * U[:, n] - 1im * k * beta * (b[1, 1] * f1 + b[1, 2] * f2 + b[1, 3] * f3 + b[1, 4] * f4)
            F2 = a[2] * M * U[:, n] - 1im * k * beta * (b[2, 1] * f1 + b[2, 2] * f2 + b[2, 3] * f3 + b[2, 4] * f4)


            # solve linear equations
            xi_new_1 = LU_Matrix_1 \ F1
            xi_new_2 = LU_Matrix_2 \ F2

            # update iteration paramter
            dU_1 = xi_old_1 - xi_new_1
            dU_2 = xi_old_2 - xi_new_2
            delta = sqrt(real(dot(dU_1, M * dU_1) + dot(dU_2, M * dU_2)))
            copyto!(xi_old_1, xi_new_1)
            copyto!(xi_old_2, xi_new_2)

            l = l + 1

            print("->")
            print(l)

        end

        # error message and break if iteration didnt converge
        if l == l_max
            println("no convergence")
            break
        end

        if isnan(delta)
            println("no convergence")
            break
        end

        # set result and calculate U(n+1)
        xi_transform_1 = Si[1, 1] * xi_old_1 + Si[1, 2] * xi_old_2
        xi_transform_2 = Si[2, 1] * xi_old_1 + Si[2, 2] * xi_old_2

        U[:, n+1] .= li_hat[1](1) * U[:, n] + li_hat[2](1) * xi_transform_1 + li_hat[3](1) * xi_transform_2

        # intermediate output
        println(" ")
        println(n)
    end

    return U

end