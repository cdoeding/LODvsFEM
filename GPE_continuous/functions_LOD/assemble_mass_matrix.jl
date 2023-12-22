function assembleGlobalMassMatrix(Tri, nodes)

    N = length(nodes)
    nodes2mesh = Array{Int}(1:N)
    M = spzeros(N, N)
    M_loc = zeros(2, 2)

    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    for i in eachindex(Tri[1, :])
        simplex .= Tri[:, i]
        detBT = abs(nodes[simplex[2]] - nodes[simplex[1]])
        M_loc .= (detBT / 6) * [2 1; 1 2]

        idx .= nodes2mesh[simplex]
        for i = 1:2
            for j = 1:2
                M[idx[i], idx[j]] += (M_loc[i, j])
            end
        end
    end

    return M

end