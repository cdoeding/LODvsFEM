function assembleGlobalPotentialMatrix(V, Tri, nodes)

    N = length(nodes)
    nodes2mesh = Array{Int}(1:N)
    M = spzeros(N, N)
    M_loc = zeros(2, 2)

    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    tau, w = get_gauss_nodes(3, 0, 1)
    phi1 = x -> 1 - x
    phi2 = x -> x

    for i in eachindex(Tri[1, :])
        simplex .= Tri[:, i]
        idx .= nodes2mesh[simplex]

        xr = nodes[simplex[2]]
        xl = nodes[simplex[1]]
        detBT = abs(xr - xl)

        V1 = V(xl + detBT * tau[1])
        V2 = V(xl + detBT * tau[2])
        V3 = V(xl + detBT * tau[3])

        M_loc[1, 1] = detBT * (w[1] * V1 * phi1(tau[1])^2 + w[2] * V2 * phi1(tau[2])^2 + w[3] * V3 * phi1(tau[3])^2)
        M_loc[2, 2] = detBT * (w[1] * V1 * phi2(tau[1])^2 + w[2] * V2 * phi2(tau[2])^2 + w[3] * V3 * phi2(tau[3])^2)
        M_loc[1, 2] = detBT * (w[1] * V1 * phi1(tau[1]) * phi2(tau[1]) + w[2] * V2 * phi1(tau[2]) * phi2(tau[2]) + w[3] * V3 * phi1(tau[3]) * phi2(tau[3]))
        M_loc[2, 1] = M_loc[1, 2]

        for i = 1:2
            for j = 1:2
                M[idx[i], idx[j]] += (M_loc[i, j])
            end
        end
    end

    return M

end

function assembleLocalPotentialMatrix(Tri_h, nodes_h, List_sub, l)

    active_fine_triangles = List_sub[:, l]
    Tri = Tri_h[:, active_fine_triangles[1]:active_fine_triangles[length(active_fine_triangles)]]
    N = length(nodes_h)
    nodes2mesh = Array{Int}(1:N)
    M = spzeros(N, N)
    M_loc = zeros(2, 2)

    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    tau, w = get_gauss_nodes(3, 0, 1)
    phi1 = x -> 1 - x
    phi2 = x -> x

    for i in eachindex(Tri[1, :])
        simplex .= Tri[:, i]
        idx .= nodes2mesh[simplex]

        xr = nodes_h[simplex[2]]
        xl = nodes_h[simplex[1]]
        detBT = abs(xr - xl)

        V1 = V(xl + detBT * tau[1])
        V2 = V(xl + detBT * tau[2])
        V3 = V(xl + detBT * tau[3])

        M_loc[1, 1] = detBT * (w[1] * V1 * phi1(tau[1])^2 + w[2] * V2 * phi1(tau[2])^2 + w[3] * V3 * phi1(tau[3])^2)
        M_loc[2, 2] = detBT * (w[1] * V1 * phi2(tau[1])^2 + w[2] * V2 * phi2(tau[2])^2 + w[3] * V3 * phi2(tau[3])^2)
        M_loc[1, 2] = detBT * (w[1] * V1 * phi1(tau[1]) * phi2(tau[1]) + w[2] * V2 * phi1(tau[2]) * phi2(tau[2]) + w[3] * V3 * phi1(tau[3]) * phi2(tau[3]))
        M_loc[2, 1] = M_loc[1, 2]

        for i = 1:2
            for j = 1:2
                M[idx[i], idx[j]] += (M_loc[i, j])
            end
        end
    end

    return M

end