function getCorrectorMatrix(T_H, p_H, T_h, p_h, Layer, List_sub, P1, S_h, M_h, B_H, B_h, alpha)

    N_h = size(p_h, 1)
    N_H = size(p_H, 1)
    NT_H = size(T_H, 2)

    C_h = P1 * M_h
    Q = spzeros(N_H, N_h)

    for l = 1:NT_H
        Rl_H, Rl_h = getRestriction(T_H, p_H, T_h, p_h, l, Layer, List_sub, B_H, B_h)
        T = getNode2MeshMatrix(T_H, p_H, l)
        A_loc = assembleLocalStiffnessMatrix(T_h, p_h, List_sub, l)
        M_V_loc = assembleLocalPotentialMatrix(T_h, p_h, List_sub, l)

        Nl_h = size(Rl_h, 1)
        Nl_H = size(Rl_H, 1)

        A = Rl_h * S_h * Rl_h'
        C = Array(Rl_H * C_h * Rl_h')

        rhs = -T * B_H * P1 * (imag(alpha) * A_loc + M_V_loc) * Rl_h'

        Y = A \ C'
        S_inv = (C * Y) \ Matrix(1.0I, Nl_H, Nl_H)


        w = zeros(2, Nl_h)
        for i = 1:2
            q = A \ Array(rhs[i, :])
            lambda = S_inv * C * q
            w[i, :] = (q - Y * lambda)'
        end

        Q += T' * sparse(w) * Rl_h

        if mod(l, 10) == 0
            println(string(l) * " patches of " * string(NT_H) * " computed")
        end

    end

    return Q

end

function getRestriction(T_H, p_H, T_h, p_h, l, Layer, List_sub, B_H, B_h)

    N_h = length(p_h)
    N_H = length(p_H)

    # coarse restriction
    active_triangles = T_H[:, findnz(Layer[:, l])[2]]
    active_nodes = Array{Bool}(zeros(N_H, 1))

    for i = 1:2
        for j in eachindex(active_triangles[1, :])
            active_nodes[active_triangles[i, j]] = true
        end
    end

    active_nodes = B_H * active_nodes
    N_H_active = Int(sum(active_nodes))
    Rl_H = spzeros(N_H_active, N_H)

    index = 1
    for i = 1:N_H
        if active_nodes[i] == true
            Rl_H[index, i] = 1
            index = index + 1
        end
    end

    # # fine restriction 
    active_fine_layer_init = reshape(List_sub[:, findnz(Layer[:, l])[2]], 1, :)
    active_fine_layer = Vector{Int64}(undef, size(active_fine_layer_init, 2))
    copyto!(active_fine_layer, active_fine_layer_init)
    active_fine_triangles = T_h[:, active_fine_layer[1]:active_fine_layer[length(active_fine_layer)]]
    active_fine_nodes = Array{Bool}(zeros(N_h, 1))

    for i = 1:2
        for j in eachindex(active_fine_triangles[1, :])
            active_fine_nodes[active_fine_triangles[i, j]] = true
        end
    end

    active_fine_nodes = B_h * active_fine_nodes
    N_h_active = Int(sum(active_fine_nodes))
    Rl_h = spzeros(N_h_active, N_h)

    index = 1
    for i = 1:N_h
        if active_fine_nodes[i] == true
            Rl_h[index, i] = 1
            index = index + 1
        end
    end

    return Rl_H, Rl_h

end

function getNode2MeshMatrix(T_H, p_H, l)

    N_H = size(p_H, 1)
    T = spzeros(2, N_H)

    for k = 1:2
        T[k, t_H[k, l]] = 1
    end

    return T

end




























