function assemble_tensor(W_Iptr, W_J, W_K, W_val, VLOD, N_H, t_h, p_h)

    nodes2mesh_h = zeros(Int64, length(p_h), 1)
    nodes2mesh_h[2:length(p_h)-1] = 1:length(p_h)-2
    N_h = size(t_h, 2)

    max_alloc = 1
    for index = 1:N_H

        start = W_Iptr[index]
        finish = W_Iptr[index+1] - 1

        max_loc = finish - start + 1
        if (max_loc > max_alloc)
            max_alloc = max_loc
        end
    end

    tau = [0.5 (sqrt(3 / 5) + 1.0) / 2.0 (1.0 - sqrt(3 / 5)) / 2.0]
    w = [4 / 9
        5 / 18
        5 / 18]

    max_number = 120
    quadrature_points = rand(max_number, 3)
    active_nodes = Vector{Int64}(undef, 2000)

    phi1 = 1 .- tau
    phi = [phi1; tau]

    idx = size(VLOD, 2)
    ta = time()
    n_it = 0

    for n = 1:N_h
        n_it += 1
        if (mod(n_it, 10000) == 0)
            println(n_it, " ", time() - ta, "percentage done ", n_it / N_h)
        end

        simplex = t_h[:, n]
        meshidx = nodes2mesh_h[simplex]
        Vertices = p_h[simplex]
        J = Vertices[2] - Vertices[1]

        idx = 1
        for i = 1:2
            if (meshidx[i] != 0)
                active_nodes[idx:idx+length(VLOD[:, meshidx[i]].nzind)-1] = VLOD[:, meshidx[i]].nzind
                idx += length(VLOD[:, meshidx[i]].nzind)
            end
        end

        a_n = active_nodes[1:idx-1]
        sort!(a_n)
        unique!(a_n)

        if (!isempty(a_n))
            len = length(a_n)
            DetJ = (det(J))

            for i = 1:2
                if (meshidx[i] == 0)
                    meshidx[i] = idx
                end
            end

            quadrature_points[1:len, :] .= VLOD[a_n, meshidx] * phi
            quadrature_rule(quadrature_points, w, len, DetJ, a_n, W_Iptr, W_J, W_K, W_val, max_alloc)
        end
    end
end

function quadrature_rule(quadrature_points, w, len, DetJ, a_n, W_Iptr, W_J, W_K, W_val, max_alloc)
    Possible_IdxJ = zeros(Int64, max_alloc)

    for i = 1:len
        i_glob = a_n[i]

        start = W_Iptr[i_glob]
        finish = W_Iptr[i_glob+1] - 1

        Possible_IdxJ[1:(finish-start+1)] = W_J[start:finish]
        Possible_IdxJ[finish-start+2:end] .= 2^100

        f1 = quadrature_points[i, 1] * w[1] * DetJ
        f2 = quadrature_points[i, 2] * w[2] * DetJ
        f3 = quadrature_points[i, 3] * w[3] * DetJ

        Possible_IdxJ[1:(finish-start+1)] = W_J[start:finish]
        Possible_IdxJ[finish-start+2:end] .= 2^100

        for j = i:len
            j_glob = a_n[j]
            sub_j = searchsorted(Possible_IdxJ, j_glob)

            g1 = f1 * quadrature_points[j, 1]
            g2 = f2 * quadrature_points[j, 2]
            g3 = f3 * quadrature_points[j, 3]
            if (!isempty(sub_j))

                Possible_IdxK = W_K[start.+sub_j.-1]
                for k = j:len
                    k_glob = a_n[k]
                    sub_k = searchsorted(Possible_IdxK, k_glob)

                    val = g1 * quadrature_points[k, 1] + g2 * quadrature_points[k, 2] + g3 * quadrature_points[k, 3]
                    W_val[(start+sub_j[1]-2).+sub_k] .+= val
                end
            end
        end
    end

end
