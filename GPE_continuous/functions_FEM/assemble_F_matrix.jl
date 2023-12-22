function assemble_F_matrix(v, Tri, nodes, nodes2mesh)

    N = length(nodes) - 1
    F = spzeros(N - 1, N - 1) * 0im
    F_loc = zeros(2, 2) * 0im

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

        if idx[1] == 0
            u = x -> v[idx[2]] * phi2(x)
        elseif idx[2] == 0
            u = x -> v[idx[1]] * phi1(x)
        else
            u = x -> v[idx[1]] * phi1(x) + v[idx[2]] * phi2(x)
        end

        f1 = abs(u(tau[1]))^2 * u(tau[1])
        f2 = abs(u(tau[2]))^2 * u(tau[2])
        f3 = abs(u(tau[3]))^2 * u(tau[3])

        F_loc[1, 1] = detBT * (w[1] * f1 * phi1(tau[1])^2 + w[2] * f2 * phi1(tau[2])^2 + w[3] * f3 * phi1(tau[3])^2)
        F_loc[2, 2] = detBT * (w[1] * f1 * phi2(tau[1])^2 + w[2] * f2 * phi2(tau[2])^2 + w[3] * f3 * phi2(tau[3])^2)
        F_loc[1, 2] = detBT * (w[1] * f1 * phi1(tau[1]) * phi2(tau[1]) + w[2] * f2 * phi1(tau[2]) * phi2(tau[2]) + w[3] * f3 * phi1(tau[3]) * phi2(tau[3]))
        F_loc[2, 1] = F_loc[1, 2]

        for i = 1:2
            for j = 1:2
                if ((idx[i] != 0) & (idx[j] != 0))
                    F[idx[i], idx[j]] += (F_loc[i, j])
                end
            end
        end
    end

    return M

end

function assemble_F_P2(v, Tri, nodes, nodes2mesh)

    f = zeros(size(v)) * 0im

    idx = Array{Int}(zeros(3))
    simplex = Array{Int}(zeros(3))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        idx .= nodes2mesh[simplex]

        xr = nodes[simplex[3]]
        xl = nodes[simplex[1]]
        detBT = abs(xr - xl)

        tau, w = get_gauss_nodes(6, 0, 1)
        phi1 = x -> 2 * (x - 0.5) * (x - 1)
        phi2 = x -> -4 * x * (x - 1)
        phi3 = x -> 2 * x * (x - 0.5)

        if idx[1] == 0
            u = x -> v[idx[2]] * phi2(x) + v[idx[3]] * phi3(x)
        elseif idx[3] == 0
            u = x -> v[idx[1]] * phi1(x) + v[idx[2]] * phi2(x)
        else
            u = x -> v[idx[1]] * phi1(x) + v[idx[2]] * phi2(x) + v[idx[3]] * phi3(x)
        end

        f1 = abs(u(tau[1]))^2 * u(tau[1])
        f2 = abs(u(tau[2]))^2 * u(tau[2])
        f3 = abs(u(tau[3]))^2 * u(tau[3])
        f4 = abs(u(tau[4]))^2 * u(tau[4])
        f5 = abs(u(tau[5]))^2 * u(tau[5])
        f6 = abs(u(tau[6]))^2 * u(tau[6])

        if idx[1] != 0
            f[idx[1]] += detBT * (w[1] * f1 * phi1(tau[1]) + w[2] * f2 * phi1(tau[2]) + w[3] * f3 * phi1(tau[3]) + w[4] * f4 * phi1(tau[4]) + w[5] * f5 * phi1(tau[5]) + w[6] * f6 * phi1(tau[6]))
        end

        if idx[2] != 0
            f[idx[2]] += detBT * (w[1] * f1 * phi2(tau[1]) + w[2] * f2 * phi2(tau[2]) + w[3] * f3 * phi2(tau[3]) + w[4] * f4 * phi2(tau[4]) + w[5] * f5 * phi2(tau[5]) + w[6] * f6 * phi2(tau[6]))
        end

        if idx[3] != 0
            f[idx[3]] += detBT * (w[1] * f1 * phi3(tau[1]) + w[2] * f2 * phi3(tau[2]) + w[3] * f3 * phi3(tau[3]) + w[4] * f4 * phi3(tau[4]) + w[5] * f5 * phi3(tau[5]) + w[6] * f6 * phi3(tau[6]))
        end

    end

    return f

end

function assemble_F_P3(v, Tri, nodes, nodes2mesh)

    f = zeros(size(v)) * 0im

    idx = Array{Int}(zeros(4))
    simplex = Array{Int}(zeros(4))

    for i in eachindex(Tri[:, 1])
        simplex .= Tri[i, :]
        idx .= nodes2mesh[simplex]

        xr = nodes[simplex[4]]
        xl = nodes[simplex[1]]
        detBT = abs(xr - xl)

        tau, w = get_gauss_nodes(6, 0, 1)
        phi1 = x -> -9 / 2 * (x - 1 / 3) * (x - 2 / 3) * (x - 1)
        phi2 = x -> 27 / 2 * x * (x - 2 / 3) * (x - 1)
        phi3 = x -> -27 / 2 * x * (x - 1 / 3) * (x - 1)
        phi4 = x -> 9 / 2 * x * (x - 1 / 3) * (x - 2 / 3)

        if idx[1] == 0
            u = x -> v[idx[2]] * phi2(x) + v[idx[3]] * phi3(x) + v[idx[4]] * phi4(x)
        elseif idx[4] == 0
            u = x -> v[idx[1]] * phi1(x) + v[idx[2]] * phi2(x) + v[idx[3]] * phi3(x)
        else
            u = x -> v[idx[1]] * phi1(x) + v[idx[2]] * phi2(x) + v[idx[3]] * phi3(x) + v[idx[4]] * phi4(x)
        end

        f1 = abs(u(tau[1]))^2 * u(tau[1])
        f2 = abs(u(tau[2]))^2 * u(tau[2])
        f3 = abs(u(tau[3]))^2 * u(tau[3])
        f4 = abs(u(tau[4]))^2 * u(tau[4])
        f5 = abs(u(tau[5]))^2 * u(tau[5])
        f6 = abs(u(tau[6]))^2 * u(tau[6])

        if idx[1] != 0
            f[idx[1]] += detBT * (w[1] * f1 * phi1(tau[1]) + w[2] * f2 * phi1(tau[2]) + w[3] * f3 * phi1(tau[3]) + w[4] * f4 * phi1(tau[4]) + w[5] * f5 * phi1(tau[5]) + w[6] * f6 * phi1(tau[6]))
        end

        if idx[2] != 0
            f[idx[2]] += detBT * (w[1] * f1 * phi2(tau[1]) + w[2] * f2 * phi2(tau[2]) + w[3] * f3 * phi2(tau[3]) + w[4] * f4 * phi2(tau[4]) + w[5] * f5 * phi2(tau[5]) + w[6] * f6 * phi2(tau[6]))
        end

        if idx[3] != 0
            f[idx[3]] += detBT * (w[1] * f1 * phi3(tau[1]) + w[2] * f2 * phi3(tau[2]) + w[3] * f3 * phi3(tau[3]) + w[4] * f4 * phi3(tau[4]) + w[5] * f5 * phi3(tau[5]) + w[6] * f6 * phi3(tau[6]))
        end

        if idx[4] != 0
            f[idx[4]] += detBT * (w[1] * f1 * phi4(tau[1]) + w[2] * f2 * phi4(tau[2]) + w[3] * f3 * phi4(tau[3]) + w[4] * f4 * phi4(tau[4]) + w[5] * f5 * phi4(tau[5]) + w[6] * f6 * phi4(tau[6]))
        end

    end

    return f

end