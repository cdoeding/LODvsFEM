function evaluate_initial_value(u0, t_h, p_h, M_LOD, P1, Q, B_H, B_h)

    
    r = assembleGlobalRHS(u0, t_h, p_h)
    u = M_LOD' \ (B_H * (P1 + Q) * B_h' * B_h * r)

    return u

end

function assembleGlobalRHS(f, t_h, p_h)

    F = zeros(length(p_h), 1)

    tau, w = get_gauss_nodes(3, 0, 1)
    phi1 = x -> 1 - x
    phi2 = x -> x

    for i in eachindex(t_h[1, :])
        nodes = p_h[t_h[:, i]]
        detJ = nodes[2] - nodes[1]
        F[i] += w[1] * detJ * f(nodes[1] + detJ * tau[1]) * phi1(tau[1])
        +w[2] * detJ * f(nodes[1] + detJ * tau[2]) * phi1(tau[2])
        +w[3] * detJ * f(nodes[1] + detJ * tau[3]) * phi1(tau[3])

        F[i+1] += w[1] * detJ * f(nodes[1] + detJ * tau[1]) * phi2(tau[1])
        +w[2] * detJ * f(nodes[1] + detJ * tau[2]) * phi2(tau[2])
        +w[3] * detJ * f(nodes[1] + detJ * tau[3]) * phi2(tau[3])
    end

    return F

end

