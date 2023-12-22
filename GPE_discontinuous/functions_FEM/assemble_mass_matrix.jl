function assemble_mass_matrix(Tri, nodes, nodes2mesh)

    N = length(nodes) - 1
    M = spzeros(N - 1, N - 1)
    M_loc = zeros(2, 2)

    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        detBT = abs(nodes[simplex[2]] - nodes[simplex[1]])
        M_loc .= (detBT / 6) * [2 1; 1 2]

        idx .= nodes2mesh[simplex]
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

function assemble_mass_matrix_P2(Tri, nodes, nodes2mesh)

    N = size(Tri, 1)
    M = spzeros(2 * N - 1, 2 * N - 1)
    M_loc = zeros(3, 3)

    idx = Array{Int}(zeros(3))
    simplex = Array{Int}(zeros(3))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        detBT = abs(nodes[simplex[3]] - nodes[simplex[1]])
        M_loc .= (detBT / 30) * [4 2 -1; 2 16 2; -1 2 4]

        idx .= nodes2mesh[simplex]
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

function assemble_mass_matrix_P3(Tri, nodes, nodes2mesh)

    N = size(Tri, 1)
    M = spzeros(3 * N - 1, 3 * N - 1)
    M_loc = zeros(4, 4)

    idx = Array{Int}(zeros(4))
    simplex = Array{Int}(zeros(4))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        detBT = abs(nodes[simplex[4]] - nodes[simplex[1]])
        M_loc .= detBT * [8/105 33/560 -3/140 19/1680;
            33/560 27/70 -27/560 -3/140;
            -3/140 -27/560 27/70 33/560;
            19/1680 -3/140 33/560 8/105]

        idx .= nodes2mesh[simplex]
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




