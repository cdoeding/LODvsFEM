function assembleGlobalStiffnessMatrix(Tri, nodes)

    N = length(nodes)
    nodes2mesh = Array{Int}(1:N)
    A = spzeros(N, N)
    A_loc = zeros(2, 2)

    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    for i in eachindex(Tri[1, :])
        simplex .= Tri[:, i]
        detBT = abs(nodes[simplex[2]] - nodes[simplex[1]])
        A_loc .= (1 / detBT) * [1 -1; -1 1]

        idx .= nodes2mesh[simplex]
        for k = 1:2
            for j = 1:2
                A[idx[k], idx[j]] += (A_loc[k, j])
            end
        end
    end

    return A

end

function assembleLocalStiffnessMatrix(Tri_h, nodes_h, List_sub, l)

    active_fine_triangles = List_sub[:, l]
    Tri = Tri_h[:, active_fine_triangles[1]:active_fine_triangles[length(active_fine_triangles)]]
    N = length(nodes_h)
    nodes2mesh = Array{Int}(1:N)
    A = spzeros(N, N)
    A_loc = zeros(2, 2)

    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    for i in eachindex(Tri[1, :])
        simplex .= Tri[:, i]
        detBT = abs(nodes_h[simplex[2]] - nodes_h[simplex[1]])
        A_loc .= (1 / detBT) * [1 -1; -1 1]

        idx .= nodes2mesh[simplex]
        for i = 1:2
            for j = 1:2
                A[idx[i], idx[j]] += (A_loc[i, j])
            end
        end
    end

    return A

end