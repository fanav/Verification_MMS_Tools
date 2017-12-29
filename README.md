## Verification MMS Tools: WallFlow_Euler/NS/RANS-SA
This branch provides tools for the verification of Euler/NS/RANS-SA solvers via the method of manufactured solutions (MMS).

---------------------------------------------
#### -> Please Cite us upon utilization <-
The reference describing the test cases is: Coming soon     

""

Farshad Navah (farshad.navah@mail.mcgill.ca) ; Siva Nadarajah 

McGill University

---------------------------------------------
### Description of cases:

- MS-1: Wall-bounded inviscid flow on curved domain 
- MS-2: Wall-bounded laminar flow on curved domain 
- MS-3: Realistic boundary layer RANS-SA [1] flow on Cartesian domain, Eça et al. (2007) [2]
- MS-4: Realistic boundary layer RANS-SA [1] flow on Cartesian domain, Oliver et al. (2012) [3]

[1] Allmaras, S. R., Johnson, F. T., Spalart, P. R., "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model", ICCFD7-1902.

[2] Eça, L. , Hoekstra, M., Hay, A., Pelletier, D., "Verification of RANS solvers with manufactured solutions", doi:10.1007/s00366-007-0067-9.

[3] Oliver, T., Estacio-Hiroms, K., Malaya, N., Carey, G., "Manufactured Solutions for the Favre-Averaged Navier-Stokes Equations with Eddy-Viscosity Turbulence Models", doi:10.2514/6.2012-80, AIAA 2012-0080.

### Description of files:

- IPython3 notebook: allows to generate the manufactured solution fields as well as the forcing functions for each of the cases. For the above MSs, the resulting expressions are already injected in the .c file.

- C file: provides a structure for the implementation of the verification methodology of MMS in any (including high-order) CFD code. One call to the routine with proper identification tags for a given solution variable and one of its fields (U, dUdx, dUdy, etc. or FFU (forcing function of x-momentum)) allows to initiate a series of recursive calls to the same routine (and its sub-routines) until all the terms of the demanded field are computed and the answer is reported back to the calling instance. This allows for a facilitated and point-wise application of the verification methodology to the solver. The calls are necessary for balancing the residuals, computing the boundary conditions and initializing the solution only. Furthermore, the structure of the routines is designed to allow an optimal debugging experience. Please read the comment section in the "MS_qnt_2D" routine for more details.
