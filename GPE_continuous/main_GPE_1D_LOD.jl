# solve the 1D Gross-Piteaevskii equation
# using the space-time finite element method:
# cG(2)-method for time discretization
# LOD for spatial discretization
#
# iu_t = -u_xx + Vu + |u^2| u on [x_a,x_b] \times [0,T], u(0) = u0

using SparseArrays
using LinearAlgebra
using Preconditioners
using IterativeSolvers
using MAT

include("./functions_LOD/get_triangulation_LOD.jl")
include("./functions_LOD/get_boundary_nodes.jl")
include("./functions_LOD/assemble_stiffness_matrix.jl")
include("./functions_LOD/assemble_mass_matrix.jl")
include("./functions_LOD/assemble_potential_matrix.jl")
include("./functions_LOD/assemble_corrector.jl")
include("./functions_LOD/preallocate_tensor.jl")
include("./functions_LOD/assemble_tensor.jl")
include("./functions_LOD/get_gauss_nodes.jl")
include("./functions_LOD/evaluate_initial_value.jl")
include("./functions_LOD/assemble_nonlinear_term_by_tensor.jl")
include("./functions_LOD/solve_GPE_1D_q2.jl")

#----- system parameter -----#
alpha = 1im; # diffusion (1im for Schroedinger)
beta = 100; # repulsing parameter
V = x -> 10 * x .^ 2 # trapping potential
T = 0.4; # end time
x_a = -15; # left boundary
x_b = 15;  # right boundary

#----- numerical parameter -----#
H_level = 7;
NH = 2^H_level; # number of coarse elements 
Nh = 2^16; # number of fine elements
l = H_level + 5; # number of layers (oversampling)

Nt = 200; # number of time steps
q = 2; # degree in time
threshold = 10^(-10); # threshold for fixed point iteration
l_max = 200; # maximum iterations
k = T / Nt;

#----- solve -----#
# generate meshes for LOD
println("start triangulation")
t_H, p_H, t_h, p_h, Actual_N, Layer, List_sub, P1, P0 = get_triangulation_LOD(NH, Nh, l, x_a, x_b);

println("start boundary restriction")
B_H = getBoundaryNodes(p_H, x_a, x_b);
B_h = getBoundaryNodes(p_h, x_a, x_b);
Bd_H = getBoundaryRestriction(B_H);
Bd_h = getBoundaryRestriction(B_h);

# get standard system matricies (on fine grid)
println("start system matrix")
tid = time();
A_h = assembleGlobalStiffnessMatrix(t_h, p_h);
M_h = assembleGlobalMassMatrix(t_h, p_h);
M_V = assembleGlobalPotentialMatrix(V, t_h, p_h);

# compute corrector and LOD matricies
Q = getCorrectorMatrix(t_H, p_H, t_h, p_h, Layer, List_sub, P1, imag(alpha) * A_h + M_V, M_h, B_H, B_h, alpha);

A_LOD = Bd_H * (P1 + Q) * A_h * (P1 + Q)' * Bd_H';
M_LOD = Bd_H * (P1 + Q) * M_h * (P1 + Q)' * Bd_H';
M_V_LOD = Bd_H * (P1 + Q) * M_V * (P1 + Q)' * Bd_H';

# LOD tranformation matrix with boundary restriction
V_LOD = Bd_H * (P1 + Q) * Bd_h'

# compute tensor
println("start tensor")
W_Iptr, W_J, W_K, W_val = preallocate_tensor(V_LOD, M_LOD);
assemble_tensor(W_Iptr, W_J, W_K, W_val, V_LOD, length(p_H) - 2, t_h, p_h)
println("done")

file = matopen("./ground_state.mat")
u0_int = read(file, "u")
close(file)
u0_h = [0; u0_int; 0]
u0_vec = zeros(size(A_LOD, 1));
cg!(u0_vec, M_LOD, V_LOD * (Bd_h * M_h * Bd_h' * Bd_h * u0_h), reltol=10^-11);
time_offline = time() - tid;

# time stepping
println("start time stepping")
tid = time();
u = solve_GPE_1D_q2(u0_vec, k, Nt, A_LOD, M_LOD, M_V_LOD, W_Iptr, W_J, W_K, W_val, alpha, beta, threshold, l_max);
time_online = time() - tid;
println("done")
println(time_online)

# #----- postprocessing -----#
u_sol_h = (P1 + Q)' * Bd_H' * u;
u_h = u_sol_h[:, Nt+1];

# reference solution
file = matopen("./reference_solution.mat");
u_ref_sol = read(file, "u");
u_ref = [0; u_ref_sol[:, Nt+1]; 0];
close(file);

# compute error on fine mesh
e = u_h - u_ref;
err_H1 = abs(sqrt(e' * (A_h + M_h) * e));
err_L2 = abs(sqrt(e' * M_h * e));

dofs = size(A_LOD, 1);
H = (x_b - x_a) / NH;
h = (x_b - x_a) / Nh;

#----- save to MAT file-----#
file_name = "LOD_Hlevel" * string(H_level);
file = matopen(file_name * ".mat", "w")
write(file, "u_h", u_h)
write(file, "u_sol_h", u_sol_h)
write(file, "e", e)
write(file, "err_H1", err_H1)
write(file, "err_L2", err_L2)
write(file, "dofs", dofs)
write(file, "time_online", time_online)
write(file, "time_offline", time_offline)
write(file, "H", H)
write(file, "h", h)
write(file, "NH", NH)
write(file, "Nh", Nh)
write(file, "l", l)
write(file, "T", T)
write(file, "Nt", Nt)
write(file, "threshold", threshold)
write(file, "q", q)
close(file)
