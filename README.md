# LODvsFEM
Code used in the paper "Localized orthogonal decomposition methods vs. classical FEM for the Gross-Pitaevskii equation"

for computations with continuous potential V_1 go to /GPE_continuous

for computations with discontinuous potential V_2 go to /GPE_discontinuous

then execute the following commands e.g. in a Julia REPL:

1. compute initial value
   
	include("main_compute_groundstate.jl")
	
3. compute reference solution
   
	include("main_compute_reference_solution.jl")
	
5. cG-P1 solution (set desired "H_level" value for step size H)
   
	include("main_GPE_1D_P1_FEM.jl")
	
7. cG-P2 solution (set desired "H_level" value for step size H)
   
	include("main_GPE_1D_P2_FEM.jl")
	
9. cG-P3 solution (set desired "H_level" value for step size H)
    
	include("main_GPE_1D_P3_FEM.jl")
	
11. cG-LOD solution (set desired "H_level" value for step size H)
    
	include("main_GPE_1D_LOD.jl")
