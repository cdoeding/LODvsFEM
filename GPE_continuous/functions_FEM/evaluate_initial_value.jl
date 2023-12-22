function evaluate_initial_value(NH, Nh)

    file = matopen("./ground_state.mat")
    u0 = read(file, "u")
    close(file)

    ratio = Int(Nh / NH)
    u0_vec = u0[ratio:ratio:length(u0)-ratio+1]

    return u0_vec
end

function evaluate_initial_value_P2(NH, Nh)

    file = matopen("./ground_state.mat")
    u0 = read(file, "u")
    close(file)

    ratio = Int(Nh / (2 * NH))
    u0_vec = u0[ratio:ratio:length(u0)-ratio+1]

    return u0_vec
end

function evaluate_initial_value_P3(Tri_H, nodes_H, NH, Nh)

    file = matopen("./ground_state.mat")
    u0_int = read(file, "u")
    close(file)

    u0 = [0; u0_int; 0]
    u0_vec = zeros(length(nodes_H), 1)
    x_a = nodes_H[1]
    x_b = nodes_H[end]

    Tri_h, nodes_h, nodes2mesh_h = get_triangulation(Nh, x_a, x_b)

    ratio = Int(Nh / NH)
    i1 = Int(floor(ratio / 3))
    i2 = Int(floor(2 * ratio / 3))

    for k in eachindex(Tri_H[:, 1])

        # left interval point
        u0_vec[3*k-2] = u0[(k-1)*ratio+1]

        # left inner point
        index1 = (k - 1) * ratio + 1 + i1
        index2 = (k - 1) * ratio + 1 + i1 + 1

        xl = nodes_h[index1]
        xr = nodes_h[index2]
        phi1 = x -> (xr - x) / (xr - xl)
        phi2 = x -> (x - xl) / (xr - xl)
        u0_vec[3*k-1] = u0[index1] * phi1(nodes_H[3*k-1]) + u0[index2] * phi2(nodes_H[3*k-1])

        # right inner point
        index1 = (k - 1) * ratio + 1 + i2
        index2 = (k - 1) * ratio + 1 + i2 + 1

        xl = nodes_h[index1]
        xr = nodes_h[index2]
        phi1 = x -> (xr - x) / (xr - xl)
        phi2 = x -> (x - xl) / (xr - xl)
        u0_vec[3*k] = u0[index1] * phi1(nodes_H[3*k]) + u0[index2] * phi2(nodes_H[3*k])

    end

    return u0_vec[2:end-1]
end

function evaluate_initial_value_for_GS(u0, M, Tri, nodes, nodes2mesh)

    N = length(nodes) - 1
    rhs = zeros(N - 1, 1) * im
    idx = Array{Int}(zeros(2))
    simplex = Array{Int}(zeros(2))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        idx .= nodes2mesh[simplex]
        xr = nodes[simplex[2]]
        xl = nodes[simplex[1]]
        BT = xr - xl
        detBT = abs(xr - xl)
        tau, w = get_gauss_nodes(3, 0, 1)
        phi1 = x -> 1 - x
        phi2 = x -> x

        if idx[1] != 0
            rhs[idx[1]] += detBT * (w[1] * u0(BT * tau[1] + xl) * phi1(tau[1]) + w[2] * u0(BT * tau[2] + xl) * phi1(tau[2]) + w[3] * u0(BT * tau[3] + xl) * phi1(tau[3]))
        end

        if idx[2] != 0
            rhs[idx[2]] += detBT * (w[1] * u0(BT * tau[1] + xl) * phi2(tau[1]) + w[2] * u0(BT * tau[2] + xl) * phi2(tau[2]) + w[3] * u0(BT * tau[3] + xl) * phi2(tau[3]))
        end
    end

    u0_vec = M \ rhs

    return u0_vec
end