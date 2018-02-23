## Verification MMS Tools: WallFlow_Euler/NS/RANS-SA
This branch provides tools for the verification of Euler/NS/RANS-SA solvers via the method of manufactured solutions (MMS).

---------------------------------------------
#### -> Please Cite us upon utilization <-

The references describing the test cases are: 

[arXiv:1712.09478] (arXiv_1712.09478.pdf) (https://arxiv.org/abs/1712.09478)

"A comprehensive high-order solver verification methodology for free fluid flows"

Farshad Navah (farshad.navah@mail.mcgill.ca) ; Siva Nadarajah 

McGill University


[arXiv:1801.00021] (arXiv_1801.00021.pdf) (https://arxiv.org/abs/1801.00021)

"On the verification of CFD solvers of all orders of accuracy on curved wall-bounded domains and for realistic RANS flows"

Farshad Navah (farshad.navah@mail.mcgill.ca) ; Siva Nadarajah 

McGill University

---------------------------------------------
### Description of cases:

#### FREE FLOWS
- MS-1: EULER Verification case subsonic                                       
        the same as MS-1 in [1] and in branch FreeFlow_Euler/NS/RANS-SA       
                                                                             
- MS-2: EULER Verification case supersonic on curved domain                   
        the same as MS-2 in [1] and in branch FreeFlow_Euler/NS/RANS-SA       
                                                                             
- MS-3: NS Verification case                                                  
        the same as MS-3 in [1] and in branch FreeFlow_Euler/NS/RANS-SA       
                                                                             
- MS-4: RANS-SA Verification case, original SA model [4]                         
        the same as MS-4 in [1] and in branch FreeFlow_Euler/NS/RANS-SA       
        the same as MS-1 in [2] and in branch FreeFlow_RANS-SA                
                                                                             
- MS-5: RANS-SA Verification case, negative SA model [4]                         
        the same as MS-5 in [1] and in branch FreeFlow_Euler/NS/RANS-SA       
        the same as MS-2 in [2] and in branch FreeFlow_RANS-SA       

#### WALL-BOUNDED FLOWS        
- MS-6: EULER Verification case, slip wall-bounded curved domain              
        the same as MS-1 in [3] and in branch WallFlow_Euler/NS/RANS-SA       
                                                                             
- MS-7: Navier-Stokes Verification case, no-slip wall-bounded curved domain   
        the same as MS-2 in [3] and in branch WallFlow_Euler/NS/RANS-SA       
                                                                             
- MS-8: Realistic RANS-SA Verification case, Eça et al. (2007) [5]                
        the same as MS-3 in [3] and in branch WallFlow_Euler/NS/RANS-SA       
                                                                             
- MS-9: Realistic RANS-SA Verification case, Oliver et al. (2012) [6]             
        the same as MS-4 in [3] and in branch WallFlow_Euler/NS/RANS-SA       

References:
[1] Navah, F., Nadarajah, S. "A comprehensive high-order solver verification methodology for free fluid flows", arXiv:1712.09478.

[2] Navah, F., Nadarajah, S. "Motivations and methods of verification for high-order RANS solvers and solutions", AIAA2017-3290.

[3] Navah, F., Nadarajah, S. "On the verification of CFD solvers of all orders of accuracy on curved wall-bounded domains and for realistic RANS flows", arXiv:1801.00021

[4] Allmaras, S. R., Johnson, F. T., Spalart, P. R., "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model", ICCFD7-1902.

[5] Eça, L. , Hoekstra, M., Hay, A., Pelletier, D., "Verification of RANS solvers with manufactured solutions", doi:10.1007/s00366-007-0067-9.

[6] Oliver, T., Estacio-Hiroms, K., Malaya, N., Carey, G., "Manufactured Solutions for the Favre-Averaged Navier-Stokes Equations with Eddy-Viscosity Turbulence Models", doi:10.2514/6.2012-80, AIAA 2012-0080.

### Description of files:

- IPython3 notebook: allows to generate the manufactured solution fields as well as the forcing functions for each of the cases. For the above MSs, the resulting expressions are already injected in the .c file.

- C file: provides a structure for the implementation of the verification methodology of MMS in any (including high-order) CFD code. One call to the routine with proper identification tags for a given solution variable and one of its fields (U, dUdx, dUdy, etc. or FFU (forcing function of x-momentum)) allows to initiate a series of recursive calls to the same routine (and its sub-routines) until all the terms of the demanded field are computed and the answer is reported back to the calling instance. This allows for a facilitated and point-wise application of the verification methodology to the solver. The calls are necessary for balancing the residuals, computing the boundary conditions and initializing the solution only. Furthermore, the structure of the routines is designed to allow an optimal debugging experience. Please read the comment section in the "MS_qnt_2D" routine for more details.
