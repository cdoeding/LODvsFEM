# solve the 1D Gross-Piteaevskii equation for reference
# using the space-time finite element method:
# cG(2)-method for time discretization
# P1 Lagrage FEM for spatial discretization
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
V = x -> 10 * x .^ 2 # trapping potential
T = 0.4; # end time
x_a = -15; # left boundary
x_b = 15;  # right boundary

#----- numerical parameter -----#
Nt = 200; # number of time steps
H_level = 16;
Nx = 2^H_level; # (coarse) spatial elements
q = 2; # degree in time
threshold = 10^(-10); # threshold for fixed point iteration
l_max = 200; # maximum iterations

#----- solve -----#
# read intial value from file
file = matopen("./ground_state.mat")
u0_vec = read(file, "u")
close(file)

# generate mesh
Tri, nodes, nodes2mesh = get_triangulation(Nx, x_a, x_b);

# get system matricies
println("assemble system matrix")
A = assemble_stiffness_matrix(Tri, nodes, nodes2mesh);
M = assemble_mass_matrix(Tri, nodes, nodes2mesh);
M_V = assemble_potential_matrix(V, Tri, nodes, nodes2mesh);

# time stepping
println("start time stepping")
u = solve_GPE_1D_q2(A, M, M_V, alpha, beta, V, u0_vec, T, Nt, Tri, nodes, nodes2mesh, threshold, l_max);
println("done")

#----- postprocessing -----#
dofs = size(A, 1);
H = (x_b - x_a) / Nx;

#----- save to MAT file-----#
file_name = "reference_solution";
file = matopen(file_name * ".mat", "w")
write(file, "u", u)
write(file, "dofs", dofs)
write(file, "H", H)
write(file, "NH", Nx)
write(file, "T", T)
write(file, "Nt", Nt)
write(file, "threshold", threshold)
write(file, "q", q)
close(file)