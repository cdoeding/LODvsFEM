function getBoundaryNodes(p, x_a, x_b)

    N = length(p)
    B = spzeros(N, N)

    for i = 1:N
        if (abs(p[i] - x_a) > eps(10.0) && abs(p[i] - x_b) > eps(10.0))
            B[i, i] = 1
        end
    end

    return B

end

function getBoundaryRestriction(B_H)

    N_H = size(B_H, 1)
    active_nodes = diag(B_H)
    N_H_active = Int(sum(active_nodes))
    B = spzeros(N_H_active, N_H)

    index = 1
    for i = 1:N_H
        if active_nodes[i] != 0
            B[index, i] = 1
            index += 1
        end
    end

    return B

end