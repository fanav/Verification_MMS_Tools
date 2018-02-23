## Verification MMS Tools: FreeFlow_RANS-SA
This branch provides tools for the verification of RANS-SA solvers via the method of manufactured solutions (MMS).

---------------------------------------------
#### -> Please Cite us upon utilization <-
The reference describing the test cases is: [AIAA2017-3290](AIAA2017-3290.pdf) (https://arc.aiaa.org/doi/abs/10.2514/6.2017-3290 )     

"Motivations and methods of verification for high-order RANS solvers and solutions"

Farshad Navah (farshad.navah@mail.mcgill.ca) ; Siva Nadarajah 

McGill University

---------------------------------------------
### Description of cases:

- MS-1: Free turbulent flow with RANS and the *original* Spalart-Allmaras turbulence model¹ on Cartesian domain / The same as MS-4 in the "FreeFlow_Euler/NS/RANS-SA" branch
- MS-2: Free turbulent flow with RANS and the *negative* Spalart-Allmaras turbulence model¹ on Cartesian domain / The same as MS-5 in the "FreeFlow_Euler/NS/RANS-SA" branch

¹Allmaras, S. R., Johnson, F. T., Spalart, P. R., "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model". ICCFD7-1902.

### Description of files:

- IPython3 notebook: allows to generate the manufactured solution fields as well as the forcing functions for each of the cases. For the above MSs, the resulting expressions are already injected in the .c file.

- C file: provides a structure for the implementation of the verification methodology of MMS in any (including high-order) CFD code. One call to the routine with proper identification tags for a given solution variable and one of its fields (U, dUdx, dUdy, etc. or FFU (forcing function of x-momentum)) allows to initiate a series of recursive calls to the same routine (and its sub-routines) until all the terms of the demanded field are computed and the answer is reported back to the calling instance. This allows for a facilitated and point-wise application of the verification methodology to the solver. The calls are necessary for balancing the residuals, computing the boundary conditions and initializing the solution only. Furthermore, the structure of the routines is designed to allow an optimal debugging experience. Please read the comment section in the "MS_qnt_2D" routine for more details.
