# solve the 1D Gross-Piteaevskii equation
# using the space-time finite element method:
# cG(2)-method for time discretization
# P2 Lagrage FEM for spatial discretization
#
# iu_t = -u_xx + Vu + |u^2| u on [x_a,x_b] \times [0,T], u(0) = u0

using SparseArrays
using LinearAlgebra
using MAT

include("./functions_FEM/get_triangulation.jl")
include("./functions_FEM/get_gauss_nodes.jl")
include("./functions_FEM/assemble_stiffness_matrix.jl")
include("./functions_FEM/assemble_mass_matrix.jl")
include("./functions_FEM/assemble_potential_matrix.jl")
include("./functions_FEM/evaluate_initial_value.jl")
include("./functions_FEM/assemble_F.jl")
include("./functions_FEM/solve_GPE_1D_q2.jl")
include("./functions_FEM/polongate_solution.jl")

#----- system parameter -----#
alpha = 1im; # diffusion (1im for Schroedinger)
beta = 100; # repulsing parameter
V = x -> 10 * (x < 0) * x .^ 2 + 100 * (x >= 5) # trapping potential
T = 0.4; # end time
x_a = -15; # left boundary
x_b = 15;  # right boundary

#----- numerical parameter -----#
Nt = 200; # number of time steps
H_level = 7;
Nx = 2^H_level; # (coarse) spatial elements
q = 2; # degree in time
threshold = 10^(-10); # threshold for fixed point iteration
l_max = 200; # maximum iterations

#----- solve -----#
# read intial value
Nh = 2^16 # (coarse) spatial elements for intial value
u0_vec = evaluate_initial_value_P2(Nx, Nh)

# generate mesh
Tri, nodes, nodes2mesh = get_triangulation_P2(Nx, x_a, x_b);

# get system matricies
println("assemble system matrix")
tid = time();
A = assemble_stiffness_matrix_P2(Tri, nodes, nodes2mesh);
M = assemble_mass_matrix_P2(Tri, nodes, nodes2mesh);
M_V = assemble_potential_matrix_P2(V, Tri, nodes, nodes2mesh);
time_offline = time() - tid;

# time stepping
println("start time stepping")
tid = time();
u = solve_GPE_1D_q2_P2(A, M, M_V, alpha, beta, V, u0_vec, T, Nt, Tri, nodes, nodes2mesh, threshold, l_max);
time_online = time() - tid;
println("done")
println(time_online)

#----- postprocessing -----#
Tri_h, nodes_h, nodes2mesh_h = get_triangulation(Nh, x_a, x_b);
u_h = polongate_solution_P2(u[:, Nt+1], Tri, nodes, nodes2mesh, nodes_h);

A_h = assemble_stiffness_matrix(Tri_h, nodes_h, nodes2mesh_h);
M_h = assemble_mass_matrix(Tri_h, nodes_h, nodes2mesh_h);

# reference solution
file = matopen("./reference_solution.mat")
u_ref_sol = read(file, "u")
u_ref = u_ref_sol[:, Nt+1]
close(file)

# compute error on fine mesh
e = u_h - u_ref;
err_H1 = abs(sqrt(e' * (A_h + M_h) * e));
err_L2 = abs(sqrt(e' * M_h * e));

dofs = size(A, 1);
H = (x_b - x_a) / Nx;
h = (x_b - x_a) / Nh;

#----- save to MAT file-----#
file_name = "P2_Hlevel" * string(H_level);
file = matopen(file_name * ".mat", "w")
write(file, "u_h", u_h)
write(file, "u_H", u)
write(file, "e", e)
write(file, "err_H1", err_H1)
write(file, "err_L2", err_L2)
write(file, "dofs", dofs)
write(file, "time_online", time_online)
write(file, "time_offline", time_offline)
write(file, "H", H)
write(file, "h", h)
write(file, "NH", Nx)
write(file, "Nh", Nh)
write(file, "T", T)
write(file, "Nt", Nt)
write(file, "threshold", threshold)
write(file, "q", q)
close(file)
