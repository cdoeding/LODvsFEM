# solve the 1D Gross-Piteaevskii EVP
# using the  P1 Lagrage FEM
# Sobolev gradient descent
# lambda*u = -u_xx + Vu + |u^2| u on [x_a,x_b]

using SparseArrays
using LinearAlgebra
using MAT

include("./functions_FEM/get_triangulation.jl")
include("./functions_FEM/assemble_stiffness_matrix.jl")
include("./functions_FEM/assemble_mass_matrix.jl")
include("./functions_FEM/assemble_potential_matrix.jl")
include("./functions_FEM/evaluate_initial_value.jl")
include("./functions_FEM/assemble_F_matrix.jl")
include("./functions_FEM/solve_GPE_1D_q2.jl")
include("./functions_FEM/get_gauss_nodes.jl")

#----- system parameter -----#
alpha = 1im; # diffusion (1im for Schroedinger)
beta = 100; # repulsing parameter
V = x -> x .^ 2 # trapping potential
x_a = -15; # left boundary
x_b = 15;  # right boundary
u0(x) = exp.(-x .^ 2); # intial value for gradient descent

#----- numerical parameter -----#
H_level = 16;
Nx = 2^H_level; # (coarse) spatial elements
tau = 1; # gradient flow time step
tol = 10^(-12); # threshold for fixed point iteration
l_max = 20000; # maximum iterations

#----- solve -----#
# generate mesh
Tri, nodes, nodes2mesh = get_triangulation(Nx, x_a, x_b);

# get system matricies
println("assemble system matrix")
A = assemble_stiffness_matrix(Tri, nodes, nodes2mesh);
M = assemble_mass_matrix(Tri, nodes, nodes2mesh);
M_V = assemble_potential_matrix(V, Tri, nodes, nodes2mesh);

# get initial value
global u = evaluate_initial_value_for_GS(u0, M, Tri, nodes, nodes2mesh);
global F = assemble_F_matrix(u, Tri, nodes, nodes2mesh)
global E = 0.5 * (u'*(imag(alpha)*A*u+M_V*u+(beta/2)*F*u))[1]

# gradient flow approximation
delta = tol;
counter = 0;

while delta >= tol && counter < l_max
    global E_old = E

    Mat = (imag(alpha) * A + M_V + beta * F)
    global u = Mat \ (M * u)
    global u /= sqrt(u' * M * u)
    global F = assemble_F_matrix(u, Tri, nodes, nodes2mesh)
    global E = 0.5 * (u'*(imag(alpha)*A*u+M_V*u+(beta/2)*F*u))[1]
    global delta = abs(real(E - E_old))
    global counter += 1

    println(counter)
    println(delta)
    println(E)

end

#----- save to MAT file-----#
file_name = "ground_state";
file = matopen(file_name * ".mat", "w")
write(file, "u", u)
close(file)