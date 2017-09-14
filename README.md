## Verification_MMS_Tools
This Repository provides tools for the verification of RANS-SA solvers via the method of manufactured solutions (MMS).

#### -> Please Cite us upon utilization <-
---------------------------------------------
The reference describing the test cases is: AIAA2017-3290         

"Motivations and methods of verification for high-order RANS solvers and solutions"

Farshad Navah (farshad.navah@mail.mcgill.ca) ; Sivakumaran Nadarajah 

McGill University

#### Read More: https://arc.aiaa.org/doi/abs/10.2514/6.2017-3290
---------------------------------------------

### Description of files:

- IPython3 notebook: allows to generate the manufactured solution fields as well as the forcing functions for the verification of the compressible RANS equations (2D) along with either the original (MS-1) or the modified (MS-2) portions of the Spalart-Allmaras turbulence model¹.

- C routine: provides a structure for the implementation of the verification methodology of MMS in any (including high-order) CFD code. One call to the routine with proper identification tags for a given solution variable and one of its fields (U, dUdx, dUdy, etc. or QU (forcing function)) allows to initiate a series of recursive calls to the same routine (and its sub-routines) until all the terms of the demanded field are computed and the answer is reported back to the calling instance. This allows for a facilitated and point-wise application of the verification methodology to the code. The calls are necessary for balancing the residuals and computing the boundary conditions only.

¹Allmaras, S. R., Johnson, F. T., Spalart, P. R., "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model". ICCFD7-1902.
