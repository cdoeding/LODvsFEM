function assemble_potential_matrix(V, Tri, nodes, nodes2mesh)

    N = length(nodes) - 1
    M = spzeros(N - 1, N - 1)
    M_loc = zeros(2, 2)

    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    tau, w = get_gauss_nodes(3, 0, 1)
    phi1 = x -> 1 - x
    phi2 = x -> x

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
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
                if ((idx[i] != 0) & (idx[j] != 0))
                    M[idx[i], idx[j]] += (M_loc[i, j])
                end
            end
        end
    end

    return M

end

function assemble_potential_matrix_P2(V, Tri, nodes, nodes2mesh)

    N = size(Tri, 1)
    M = spzeros(2 * N - 1, 2 * N - 1)
    M_loc = zeros(3, 3)

    idx = Array{Int}(zeros(3))
    simplex = Array{Int}(zeros(3))

    tau, w = get_gauss_nodes(6, 0, 1)
    phi1 = x -> 2 * (x - 0.5) * (x - 1)
    phi2 = x -> -4 * x * (x - 1)
    phi3 = x -> 2 * x * (x - 0.5)

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        idx .= nodes2mesh[simplex]

        xr = nodes[simplex[3]]
        xl = nodes[simplex[1]]
        detBT = abs(xr - xl)

        V1 = V(xl + detBT * tau[1])
        V2 = V(xl + detBT * tau[2])
        V3 = V(xl + detBT * tau[3])
        V4 = V(xl + detBT * tau[4])
        V5 = V(xl + detBT * tau[5])
        V6 = V(xl + detBT * tau[6])

        M_loc[1, 1] = detBT * (w[1] * V1 * phi1(tau[1])^2 + w[2] * V2 * phi1(tau[2])^2 + w[3] * V3 * phi1(tau[3])^2 + w[4] * V4 * phi1(tau[4])^2 + w[5] * V5 * phi1(tau[5])^2 + w[6] * V6 * phi1(tau[6])^2)
        M_loc[2, 2] = detBT * (w[1] * V1 * phi2(tau[1])^2 + w[2] * V2 * phi2(tau[2])^2 + w[3] * V3 * phi2(tau[3])^2 + w[4] * V4 * phi2(tau[4])^2 + w[5] * V5 * phi2(tau[5])^2 + w[6] * V6 * phi2(tau[6])^2)
        M_loc[3, 3] = detBT * (w[1] * V1 * phi3(tau[1])^2 + w[2] * V2 * phi3(tau[2])^2 + w[3] * V3 * phi3(tau[3])^2 + w[4] * V4 * phi3(tau[4])^2 + w[5] * V5 * phi3(tau[5])^2 + w[6] * V6 * phi3(tau[6])^2)

        M_loc[1, 2] = detBT * (w[1] * V1 * phi1(tau[1]) * phi2(tau[1]) + w[2] * V2 * phi1(tau[2]) * phi2(tau[2]) + w[3] * V3 * phi1(tau[3]) * phi2(tau[3]) + w[4] * V4 * phi1(tau[4]) * phi2(tau[4]) + w[5] * V5 * phi1(tau[5]) * phi2(tau[5]) + w[6] * V6 * phi1(tau[6]) * phi2(tau[6]))
        M_loc[1, 3] = detBT * (w[1] * V1 * phi1(tau[1]) * phi3(tau[1]) + w[2] * V2 * phi1(tau[2]) * phi3(tau[2]) + w[3] * V3 * phi1(tau[3]) * phi3(tau[3]) + w[4] * V4 * phi1(tau[4]) * phi3(tau[4]) + w[5] * V5 * phi1(tau[5]) * phi3(tau[5]) + w[6] * V6 * phi1(tau[6]) * phi3(tau[6]))
        M_loc[2, 3] = detBT * (w[1] * V1 * phi2(tau[1]) * phi3(tau[1]) + w[2] * V2 * phi2(tau[2]) * phi3(tau[2]) + w[3] * V3 * phi2(tau[3]) * phi3(tau[3]) + w[4] * V4 * phi2(tau[4]) * phi3(tau[4]) + w[5] * V5 * phi2(tau[5]) * phi3(tau[5]) + w[6] * V6 * phi2(tau[6]) * phi3(tau[6]))

        M_loc[2, 1] = M_loc[1, 2]
        M_loc[3, 1] = M_loc[1, 3]
        M_loc[3, 2] = M_loc[2, 3]

        for i = 1:3
            for j = 1:3
                if ((idx[i] != 0) & (idx[j] != 0))
                    M[idx[i], idx[j]] += (M_loc[i, j])
                end
            end
        end
    end

    return M

end

function assemble_potential_matrix_P3(V, Tri, nodes, nodes2mesh)

    N = size(Tri, 1)
    M = spzeros(3 * N - 1, 3 * N - 1)
    M_loc = zeros(4, 4)

    idx = Array{Int}(zeros(4))
    simplex = Array{Int}(zeros(4))

    tau, w = get_gauss_nodes(6, 0, 1)
    phi1 = x -> -9 / 2 * (x - 1 / 3) * (x - 2 / 3) * (x - 1)
    phi2 = x -> 27 / 2 * x * (x - 2 / 3) * (x - 1)
    phi3 = x -> -27 / 2 * x * (x - 1 / 3) * (x - 1)
    phi4 = x -> 9 / 2 * x * (x - 1 / 3) * (x - 2 / 3)

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        idx .= nodes2mesh[simplex]

        xr = nodes[simplex[4]]
        xl = nodes[simplex[1]]
        detBT = abs(xr - xl)

        V1 = V(xl + detBT * tau[1])
        V2 = V(xl + detBT * tau[2])
        V3 = V(xl + detBT * tau[3])
        V4 = V(xl + detBT * tau[4])
        V5 = V(xl + detBT * tau[5])
        V6 = V(xl + detBT * tau[6])

        M_loc[1, 1] = detBT * (w[1] * V1 * phi1(tau[1])^2 + w[2] * V2 * phi1(tau[2])^2 + w[3] * V3 * phi1(tau[3])^2 + w[4] * V4 * phi1(tau[4])^2 + w[5] * V5 * phi1(tau[5])^2 + w[6] * V6 * phi1(tau[6])^2)
        M_loc[2, 2] = detBT * (w[1] * V1 * phi2(tau[1])^2 + w[2] * V2 * phi2(tau[2])^2 + w[3] * V3 * phi2(tau[3])^2 + w[4] * V4 * phi2(tau[4])^2 + w[5] * V5 * phi2(tau[5])^2 + w[6] * V6 * phi2(tau[6])^2)
        M_loc[3, 3] = detBT * (w[1] * V1 * phi3(tau[1])^2 + w[2] * V2 * phi3(tau[2])^2 + w[3] * V3 * phi3(tau[3])^2 + w[4] * V4 * phi3(tau[4])^2 + w[5] * V5 * phi3(tau[5])^2 + w[6] * V6 * phi3(tau[6])^2)
        M_loc[4, 4] = detBT * (w[1] * V1 * phi4(tau[1])^2 + w[2] * V2 * phi4(tau[2])^2 + w[3] * V3 * phi4(tau[3])^2 + w[4] * V4 * phi4(tau[4])^2 + w[5] * V5 * phi4(tau[5])^2 + w[6] * V6 * phi4(tau[6])^2)
        M_loc[1, 2] = detBT * (w[1] * V1 * phi1(tau[1]) * phi2(tau[1]) + w[2] * V2 * phi1(tau[2]) * phi2(tau[2]) + w[3] * V3 * phi1(tau[3]) * phi2(tau[3]) + w[4] * V4 * phi1(tau[4]) * phi2(tau[4]) + w[5] * V5 * phi1(tau[5]) * phi2(tau[5]) + w[6] * V6 * phi1(tau[6]) * phi2(tau[6]))
        M_loc[1, 3] = detBT * (w[1] * V1 * phi1(tau[1]) * phi3(tau[1]) + w[2] * V2 * phi1(tau[2]) * phi3(tau[2]) + w[3] * V3 * phi1(tau[3]) * phi3(tau[3]) + w[4] * V4 * phi1(tau[4]) * phi3(tau[4]) + w[5] * V5 * phi1(tau[5]) * phi3(tau[5]) + w[6] * V6 * phi1(tau[6]) * phi3(tau[6]))
        M_loc[1, 4] = detBT * (w[1] * V1 * phi1(tau[1]) * phi4(tau[1]) + w[2] * V2 * phi1(tau[2]) * phi4(tau[2]) + w[3] * V3 * phi1(tau[3]) * phi4(tau[3]) + w[4] * V4 * phi1(tau[4]) * phi4(tau[4]) + w[5] * V5 * phi1(tau[5]) * phi4(tau[5]) + w[6] * V6 * phi1(tau[6]) * phi4(tau[6]))
        M_loc[2, 3] = detBT * (w[1] * V1 * phi2(tau[1]) * phi3(tau[1]) + w[2] * V2 * phi2(tau[2]) * phi3(tau[2]) + w[3] * V3 * phi2(tau[3]) * phi3(tau[3]) + w[4] * V4 * phi2(tau[4]) * phi3(tau[4]) + w[5] * V5 * phi2(tau[5]) * phi3(tau[5]) + w[6] * V6 * phi2(tau[6]) * phi3(tau[6]))
        M_loc[2, 4] = detBT * (w[1] * V1 * phi2(tau[1]) * phi4(tau[1]) + w[2] * V2 * phi2(tau[2]) * phi4(tau[2]) + w[3] * V3 * phi2(tau[3]) * phi4(tau[3]) + w[4] * V4 * phi2(tau[4]) * phi4(tau[4]) + w[5] * V5 * phi2(tau[5]) * phi4(tau[5]) + w[6] * V6 * phi2(tau[6]) * phi4(tau[6]))
        M_loc[3, 4] = detBT * (w[1] * V1 * phi3(tau[1]) * phi4(tau[1]) + w[2] * V2 * phi3(tau[2]) * phi4(tau[2]) + w[3] * V3 * phi3(tau[3]) * phi4(tau[3]) + w[4] * V4 * phi3(tau[4]) * phi4(tau[4]) + w[5] * V5 * phi3(tau[5]) * phi4(tau[5]) + w[6] * V6 * phi3(tau[6]) * phi4(tau[6]))

        M_loc[2, 1] = M_loc[1, 2]
        M_loc[3, 1] = M_loc[1, 3]
        M_loc[4, 1] = M_loc[1, 4]
        M_loc[3, 2] = M_loc[2, 3]
        M_loc[4, 2] = M_loc[2, 4]
        M_loc[4, 3] = M_loc[3, 4]

        for i = 1:4
            for j = 1:4
                if ((idx[i] != 0) & (idx[j] != 0))
                    M[idx[i], idx[j]] += (M_loc[i, j])
                end
            end
        end
    end

    return M

end


