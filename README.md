# ICLS

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A brief manual for the Interface-Correction Level Set (ICLS) code

Author: Anthony Ge
Last modified: Feb 13, 2018

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Source code (.f90 files)

 **main**: Layout of the code. Start or restart the simulation here.
 
 **param**: Frequently used generic/control/governing parameters.
 
 **common**: Commonly used variables. Be careful not to duplicate those names in your local subroutines.
 
 **init**: Initialization.
 
 **bound**: Domain/processor boundary conditions.
 
 **levelset**: The main level set module.
 
 **interface**: Curvature computation, etc.
 
 **ab2**: Explicit (Adams-Bashforth) computation of the RHS of NS, before projecting to the divergence-free space.
 
 **mom**: Momentum solve (convection, diffusion, and gravity).
 
 **hyperPDE**: A collection of finite-difference schemes for solving hyperbolic PDEs.
 
 **fillps**: RHS of the Poisson equation in the projection step (except for the pressure jump terms, which are in gfm).
 
 **gfm**: Calculation of the jump terms to be added on the RHS of the Poisson equation in the projection step, using the ghost fluid method.
 
 **solver**: Pressure solve using Gauss elimination in the Fourier space.
 
 **correc**: Explicit correction of velocity after the projection.
 
 **chkdt**: Check the time step.
 
 **chkdiv**: Check the divergence of velocity.
 
 **misc**: Miscellaneous.
 
 **output**: Old-fashioned output.
 
 **write_vtk**: .vtk output.
 
 **loadd**: Checkpointing (binary format).
 
 **stopwatch**: Timing.
 
 rib: IBM mask for artificial boundary inside of the computational domain. (NOT validated)
 
 contact_line: A moving contact line model. (NOT validated)


# Libraries
 
 **2decomp**: 2D pencil decomposition of the computational domain.
 
 **fft**: Fast Fourier Transform and the cosine transform of the Poisson equation.


# Makefile

 Machine-dependent compiling commands (choose from the available options or add your own). 


# run

 A short-cut execution script.
