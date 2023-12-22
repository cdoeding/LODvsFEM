function get_triangulation_LOD(NH, Nh, l, x_a, x_b)

    H = (x_b - x_a) / NH
    N_sub = Int(ceil(Nh / NH))

    p_H, t_H, p_h, t_h, List_sub, Layer, Actual_N = generate_mesh(x_a, x_b, H, N_sub, l)

    List_sub = Array{Int}(List_sub)
    P0 = spzeros(size(t_H, 2), size(t_h, 2))

    for i in eachindex(List_sub[1, :])
        for j in eachindex(List_sub[:, i])
            index = List_sub[j, i]
            if index != 0
                P0[i, index] = 1
            end
        end
    end

    P1 = get_polongation(t_H, p_H, t_h, p_h, List_sub)

    return t_H, p_H, t_h, p_h, Actual_N, Layer, List_sub, P1, P0
end

function generate_mesh(x_a, x_b, H, Nh, l)

    N_H = Int(ceil((x_b - x_a) / H))
    p_H = Array(LinRange(x_a, x_b, N_H + 1))
    t_H = zeros(Int64, 2, N_H)
    t_H[1, :] = 1:N_H
    t_H[2, :] = 2:N_H+1
    sub_p = LinRange(0, 1, Nh + 1)
    sub_t = zeros(Int64, 2, Nh)
    sub_t[1, :] = 1:Nh
    sub_t[2, :] = 2:Nh+1
    p_h = zeros(1, length(p_H) * Nh * 2)
    t_h = zeros(Int, 2, length(t_H) * Nh * 2)
    p_count = 0
    t_count = 0
    List_sub = zeros(2, size(t_H, 2))
    Expected_N = Int(l + 1)
    Layer = spzeros(Int, Expected_N * 2, size(t_H, 2))
    Actual_N = zeros(Int, 1, size(t_H, 2))

    for i in eachindex(t_H[1, :])

        X = sub_p
        sH = t_H[:, i]
        P2 = p_H[sH[2]]
        P1 = p_H[sH[1]]
        local_sub_p = (P2 - P1) .* X .+ P1
        p_h[p_count+1:p_count+length(sub_p)] = local_sub_p
        t_h[:, t_count+1:t_count+size(sub_t, 2)] = sub_t .+ p_count
        List_sub[1, i] = t_count + 1
        List_sub[2, i] = t_count + size(sub_t, 2)
        p_count += length(sub_p)
        t_count += size(sub_t, 2)
        elemIDs = max(i - l, 1):min(N_H, i + l)
        No_elem_in_layer = length(elemIDs)
        Actual_N[i] = No_elem_in_layer
        Layer[1:No_elem_in_layer, i] = elemIDs

    end

    p_h = p_h[:, 1:p_count]
    t_h = t_h[:, 1:t_count]
    idx_unique = unique(i -> p_h[i], eachindex(p_h))
    sorted = sort(p_h[idx_unique])
    ic = zeros(Int, length(p_h))

    for i in eachindex(p_h)
        ic[i] = searchsorted(sorted, p_h[i])[1]
    end

    t_h = ic[t_h]
    p_h = p_h[idx_unique]

    return p_H, t_H, p_h, t_h, List_sub, Layer, Actual_N

end

function get_polongation(t_H, p_H, t_h, p_h, List_sub)

    P = spzeros(size(p_H, 1), size(p_h, 1))

    for i in eachindex(t_H[1, :])
        coarse_nodes = t_H[:, i]
        coarse_points = p_H[coarse_nodes]
        fine_simplices = List_sub[1, i]:List_sub[2, i]
        J = coarse_points[2] - coarse_points[1]
        x_0 = coarse_points[1]

        for j = fine_simplices
            fine_nodes = t_h[:, j]
            points = p_h[fine_nodes]
            for k = 1:2
                X = (points[k] - x_0) / J
                Vals = [1 - X[1], X[1]]
                for l = 1:2
                    P[coarse_nodes[l], fine_nodes[k]] = Vals[l]
                end
            end
        end
    end

    return P
end
