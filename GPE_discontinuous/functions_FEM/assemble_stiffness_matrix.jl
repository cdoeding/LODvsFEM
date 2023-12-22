function assemble_stiffness_matrix(Tri, nodes, nodes2mesh)

    N = length(nodes) - 1
    A = spzeros(N - 1, N - 1)
    A_loc = zeros(2, 2)

    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        detBT = abs(nodes[simplex[2]] - nodes[simplex[1]])
        A_loc .= (1 / detBT) * [1 -1; -1 1]

        idx .= nodes2mesh[simplex]
        for i = 1:2
            for j = 1:2
                if ((idx[i] != 0) & (idx[j] != 0))
                    A[idx[i], idx[j]] += (A_loc[i, j])
                end
            end
        end
    end

    return A

end

function assemble_stiffness_matrix_P2(Tri, nodes, nodes2mesh)

    N = size(Tri, 1)
    A = spzeros(2 * N - 1, 2 * N - 1)
    A_loc = zeros(3, 3)

    idx = Array{Int}(zeros(3))
    simplex = Array{Int}(zeros(3))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        detBT = abs(nodes[simplex[3]] - nodes[simplex[1]])
        A_loc .= (1 / (3 * detBT)) * [7 -8 1; -8 16 -8; 1 -8 7]

        idx .= nodes2mesh[simplex]
        for i = 1:3
            for j = 1:3
                if ((idx[i] != 0) & (idx[j] != 0))
                    A[idx[i], idx[j]] += (A_loc[i, j])
                end
            end
        end
    end

    return A

end

function assemble_stiffness_matrix_P3(Tri, nodes, nodes2mesh)

    N = size(Tri, 1)
    A = spzeros(3 * N - 1, 3 * N - 1)
    A_loc = zeros(4, 4)

    idx = Array{Int}(zeros(4))
    simplex = Array{Int}(zeros(4))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        detBT = abs(nodes[simplex[4]] - nodes[simplex[1]])
        A_loc .= (1 / detBT) * [37/10 -189/40 27/20 -13/40;
            -189/40 54/5 -297/40 27/20;
            27/20 -297/40 54/5 -189/40;
            -13/40 27/20 -189/40 37/10]

        idx .= nodes2mesh[simplex]
        for i = 1:4
            for j = 1:4
                if ((idx[i] != 0) & (idx[j] != 0))
                    A[idx[i], idx[j]] += (A_loc[i, j])
                end
            end
        end
    end

    return A

end




