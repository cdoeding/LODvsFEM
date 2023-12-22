function get_triangulation(N, x_a, x_b)

    # generate nodes and elements
    nodes = Array(LinRange(x_a, x_b, N + 1))
    T = zeros(Int64, N, 2)
    T[:, 1] = 1:N
    T[:, 2] = 2:N+1

    # Dirichlet b.c.
    nodes2mesh = zeros(Int64, N + 1, 1)
    nodes2mesh[2:N] = 1:N-1

    return T, nodes, nodes2mesh

end

function get_triangulation_P2(N, x_a, x_b)

    # generate nodes and elements
    nodes = Array(LinRange(x_a, x_b, 2 * N + 1))
    T = zeros(Int64, N, 3)
    T[:, 1] = 1:2:2*N-1
    T[:, 2] = 2:2:2*N
    T[:, 3] = 3:2:2*N+1

    # Dirichlet b.c.
    nodes2mesh = zeros(Int64, 2 * N + 1, 1)
    nodes2mesh[2:2*N] = 1:2*N-1

    return T, nodes, nodes2mesh

end

function get_triangulation_P3(N, x_a, x_b)

    # generate nodes and elements
    nodes = Array(LinRange(x_a, x_b, 3 * N + 1))
    T = zeros(Int64, N, 4)
    T[:, 1] = 1:3:3*N-2
    T[:, 2] = 2:3:3*N-1
    T[:, 3] = 3:3:3*N
    T[:, 4] = 4:3:3*N+1

    # Dirichlet b.c.
    nodes2mesh = zeros(Int64, 3 * N + 1, 1)
    nodes2mesh[2:3*N] = 1:3*N-1

    return T, nodes, nodes2mesh

end