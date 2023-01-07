# FEM for quantum harmonic oscillator

This program simulates quantum harmonic oscillator using finite elements method, on square grid. 
As a result one can obtain eigenstates and eigenvalues for a given problem and study time evolution for non-stationary states (implemented with Crank-Nicholson scheme). Results can bee seen on Psi.png, where stationary states are depicted and Eigenstate.gif, where time evolution of superposition of first and second stationary state is shown.

**main.py** contains general workflow of the program- setup of parameters and meshing at the beggining, FEM local matrices creation and solution for generalised eigenproblem. In the end it implements time evoltion of the system, creating modulus squared of the wavefunctions as a function of time and expectation value of x-coordinate of that wavefunction.

**Potentials.py** contains implementation of potentials that are added to hamiltonian matrix in the **main.py**. For quantum harmonic oscillator the potential is parabolic, more potentials to be added. 

**Mesh.py** is a module containing functions needed to operate on FEM grid.
  * *create_global_numbers()* creates global numbers of subsequent nodes
  * *find_global_number()* returns global number for a given element number and locacal node number (1-4 for square grid)
  * *find_real_space_coordinate()* returns real (x,y) coordinates for a given global number
  
  **MES.py** is a module containing functions creating needed matrices and wavefunctions for FEM calculation
  * *MES()* creates a gaussian wavefunction spanned on FEM grid 
  * *MES_from_vector()* returns a 2D matrix from 1D matrix (of wavefunctions). 1D matrices are needed to solve disretized eigenproblem, but have to be brought back to 2D in order to properly and easily plot them
  * *f()* is a functions used for easier implementation of FEM -> (x,y) coordinates
  * *shape_function()* implements linear shape functions 
  * *local_overlap_matrix()* creates local overlap matrix that is later implemented in **main.py** in S matrix to form generalized eigenproblem
  * *local_kinetic_matrix()* creates local kinetic energy matrix that is later implemented in **main.py** in H matrix to form generalized eigenproblem
  * *x_operator()* is an x position operator in FEM grid, used in **main.py** to calculate x-coordinate expectation value
