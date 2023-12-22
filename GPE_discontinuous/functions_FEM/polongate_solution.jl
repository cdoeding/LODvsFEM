function polongate_solution(u, Tri, nodes, nodes2mesh, nodes_h)

    u_h = zeros(length(nodes_h), 1) * 0im
    H = nodes[2] - nodes[1]
    h = nodes_h[2] - nodes_h[1]
    fpc = Int(H / h)
    phi1 = x -> 1 - x
    phi2 = x -> x

    for i in eachindex(Tri[:, 1])
        # construct solution in coarse reference triangle
        simplex = Tri[i, :]
        index1 = nodes2mesh[simplex[1]]
        index2 = nodes2mesh[simplex[2]]
        if index1 == 0
            f = x -> u[index2] * phi2(x)
        elseif index2 == 0
            f = x -> u[index1] * phi1(x)
        else
            f = x -> u[index1] * phi1(x) + u[index2] * phi2(x)
        end

        # polongate / evaluate at fine nodes
        for k = 1:fpc
            u_h[fpc*(i-1)+k] = f((k - 1) / fpc)
        end
    end

    return u_h[2:length(u_h)-1]
end

function polongate_solution_P2(u, Tri, nodes, nodes2mesh, nodes_h)
    u_h = zeros(length(nodes_h), 1) * 0im
    H = nodes[3] - nodes[1]
    h = nodes_h[2] - nodes_h[1]
    fpc = Int(H / h)
    phi1 = x -> 2 * (x - 0.5) * (x - 1)
    phi2 = x -> -4 * x * (x - 1)
    phi3 = x -> 2 * x * (x - 0.5)

    for i in eachindex(Tri[:, 1])
        # construct solution in coarse reference triangle
        simplex = Tri[i, :]
        index1 = nodes2mesh[simplex[1]]
        index2 = nodes2mesh[simplex[2]]
        index3 = nodes2mesh[simplex[3]]
        if index1 == 0
            f = x -> u[index2] * phi2(x) + u[index3] * phi3(x)
        elseif index3 == 0
            f = x -> u[index1] * phi1(x) + u[index2] * phi2(x)
        else
            f = x -> u[index1] * phi1(x) + u[index2] * phi2(x) + u[index3] * phi3(x)
        end

        # polongate / evaluate at fine nodes
        for k = 1:fpc
            u_h[fpc*(i-1)+k] = f((k - 1) / fpc)
        end
    end

    return u_h[2:length(u_h)-1]
end

function polongate_solution_P3(u, Tri, nodes, nodes2mesh, nodes_h)
    u_h = zeros(length(nodes_h), 1) * 0im
    H = nodes[4] - nodes[1]
    h = nodes_h[2] - nodes_h[1]
    fpc = Int(H / h)
    phi1 = x -> -9 / 2 * (x - 1 / 3) * (x - 2 / 3) * (x - 1)
    phi2 = x -> 27 / 2 * x * (x - 2 / 3) * (x - 1)
    phi3 = x -> -27 / 2 * x * (x - 1 / 3) * (x - 1)
    phi4 = x -> 9 / 2 * x * (x - 1 / 3) * (x - 2 / 3)

    for i in eachindex(Tri[:, 1])
        # construct solution in coarse reference triangle
        simplex = Tri[i, :]
        index1 = nodes2mesh[simplex[1]]
        index2 = nodes2mesh[simplex[2]]
        index3 = nodes2mesh[simplex[3]]
        index4 = nodes2mesh[simplex[4]]
        if index1 == 0
            f = x -> u[index2] * phi2(x) + u[index3] * phi3(x) + u[index4] * phi4(x)
        elseif index4 == 0
            f = x -> u[index1] * phi1(x) + u[index2] * phi2(x) + u[index3] * phi3(x)
        else
            f = x -> u[index1] * phi1(x) + u[index2] * phi2(x) + u[index3] * phi3(x) + u[index4] * phi4(x)
        end

        # polongate / evaluate at fine nodes
        for k = 1:fpc
            u_h[fpc*(i-1)+k] = f((k - 1) / fpc)
        end
    end

    return u_h[2:length(u_h)-1]
end
