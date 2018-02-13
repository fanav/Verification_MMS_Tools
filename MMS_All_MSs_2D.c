#include<stdlib.h>            /* C standard library     */
#include<stdio.h>             /* Standard input output  */
#include<math.h>              /* Mathematical functions */
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#ifdef __cplusplus
extern “C”
#else
#endif

const double piMS = 3.1415926535897932384626433832795;

/********************************************************************************/
/*                              FREE FLOWS                                      */
/********************************************************************************/
/*                                                                              */
/* REFERENCE:                                                                   */
/* [1] https://arxiv.org/abs/1712.09478                                         */
/* [2] https://arc.aiaa.org/doi/abs/10.2514/6.2017-3290                         */
/*                                                                              */
/* SOURCE: https://github.com/fanav/Verification_MMS_Tools                      */
/*                                                                              */
/* MS-1: EULER Verification case subsonic                                       */
/*     : the same as MS-1 in [1] and in branch FreeFlow_Euler/NS/RANS-SA        */
/*                                                                              */
/* MS-2: EULER Verification case supersonic on curved domain                    */
/*       the same as MS-2 in [1] and in branch FreeFlow_Euler/NS/RANS-SA        */
/*                                                                              */
/* MS-3: NS Verification case                                                   */
/*       the same as MS-3 in [1] and in branch FreeFlow_Euler/NS/RANS-SA        */
/*                                                                              */
/* MS-4: RANS-SA Verification case, original SA model                           */
/*       the same as MS-4 in [1] and in branch FreeFlow_Euler/NS/RANS-SA        */
/*       the same as MS-1 in [2] and in branch FreeFlow_RANS-SA                 */
/*                                                                              */
/* MS-5: RANS-SA Verification case, negative SA model                           */
/*       the same as MS-5 in [1] and in branch FreeFlow_Euler/NS/RANS-SA        */
/*       the same as MS-2 in [2] and in branch FreeFlow_RANS-SA                 */
/*                                                                              */
extern double MS_1(int qnt_type, int qnt_id, double x, double y);
extern double MS_2(int qnt_type, int qnt_id, double x, double y);
extern double MS_3(int qnt_type, int qnt_id, double x, double y);
extern double MS_4(int qnt_type, int qnt_id, double x, double y);
extern double MS_5(int qnt_type, int qnt_id, double x, double y);
/*                                                                              */
/********************************************************************************/


/********************************************************************************/
/*                           WALL_BOUNDED FLOWS                                 */
/********************************************************************************/
/*                                                                              */
/* REFERENCE:                                                                   */
/* [3] https://arxiv.org/abs/1801.00021                                         */
/*                                                                              */
/* SOURCE: https://github.com/fanav/Verification_MMS_Tools                      */
/*                                                                              */
/* MS-6: EULER Verification case, slip wall-bounded curved domain               */
/*       the same as MS-1 in [3] and in branch WallFlow_Euler/NS/RANS-SA        */
/*                                                                              */
/* MS-7: Navier-Stokes Verification case, no-slip wall-bounded curved domain    */
/*       the same as MS-2 in [3] and in branch WallFlow_Euler/NS/RANS-SA        */
/*                                                                              */
/* MS-8: Realistic RANS-SA Verification case, Eça et al. (2007)                 */
/*       the same as MS-3 in [3] and in branch WallFlow_Euler/NS/RANS-SA        */
/*                                                                              */
/* MS-9: Realistic RANS-SA Verification case, Oliver et al. (2012)              */
/*       the same as MS-4 in [3] and in branch WallFlow_Euler/NS/RANS-SA        */
/*                                                                              */
extern double MS_6(int qnt_type, int qnt_id, double x, double y);
extern double MS_7(int qnt_type, int qnt_id, double x, double y);
extern double MS_8(int qnt_type, int qnt_id, double x, double y);
extern double MS_9(int qnt_type, int qnt_id, double x, double y);
/*                                                                              */
/********************************************************************************/



// Enum definitions would help to reduce risk of mistake when calling MS_qnt_2D

//- The enum definition for quant type
enum qnt_type {
qt_Q           = 0 ,
qt_dQdx        = 1 ,
qt_dQdy        = 2 ,
qt_d2Qdxx      = 3 ,
qt_d2Qdxy      = 4 ,
qt_d2Qdyy      = 5 ,

qt_FFQ         = 10,
qt_FFQ_inv     = 11,
qt_FFQ_vis     = 12,
qt_FFQ_SA_prod = 13,
qt_FFQ_SA_dest = 14,
qt_FFQ_SA_dist = 15,
qt_FFQ_SA_cons = 16,
};

//- The enum definition for quant variable ID
enum qnt_id_Q {
qi_Rho    = 0 ,
qi_U      = 1 ,
qi_V      = 2 ,
qi_E      = 3 ,
qi_Ntl    = 4 ,
qi_Prsr   = 9 ,
qi_Wall_d = 10,
};

//- The enum definition for quant EQuation ID
enum qnt_id_EQ {
qi_Cont  = 0 ,
qi_Momx  = 1 ,
qi_Momy  = 2 ,
qi_Ener  = 3 ,
qi_SA    = 4 ,
};

double MS_qnt_2D(int MS_no, int qnt_type, int qnt_id, double x, double y)
{
/*

      Computes manufactured quantities (values of the manufactured solution, its derivatives or values of the forcing function)
      
      Note:
      -----
           The manufactured quantities are required to initialize the solution, impose boundary conditions, balance the residual equations via forcing functions
           and monitoring/post-processing (error computation, etc.).

           The manufactured solution and its derivatives are needed to initialize the solution, to impose boundary conditions and to compute error.

           The forcing function is needed to balance the residual equation of each PDE with the following convention:
                               PDE:                      :  EQ(Q) = div.(F(Q) - Fv(Q, dQdx)) + S(Q, dQdx) = 0
                               Forcing function          : FFQ(x) = EQ(Q_MS(x))
                               Modified residual equation: REQ(Q) = div.(F(Q) - Fv(Q, dQdx)) + S(Q, dQdx) = FFQ(x)

                               where Q is the solution variable, Q_MS the manufactured solution, F the inviscid flux vector, 
                               Fv the viscous flux vector, S a source term in the PDE (e.g. source term of the Spalart-Allmaras model)
                               and FFQ refers to the forcing function.

           In each case, a single call to MS_qnt_2D with proper arguments suffices to obtain the required quantity.	        

      Arguments :  
           MS_no    -> ID number of the manufactured case 
           qnt_type -> quantity: the solution, its first or second partial derivatives or the forcing function --> see below for full description
           qnt_id   -> primitive variable (or conservation equation): rho (continuity), u (x-momentum), v (y-momentum), E (energy), ntl (SA model), etc. --> see below for full description
           x,y      -> space coordinates where the quantities are required.

           Argument, Ms_no :           
                Ms_no=1, Ms_no=2, etc. The reference article above provides a comprehensive description of each manufactured case.

           Argument, qnt_type :           
                qnt_type=0:  Q -> primitive variable
                qnt_type=1:  dQdx
                qnt_type=2:  dQdy
                qnt_type=3:  d2Qdxx
                qnt_type=4:  d2Qdxy
                qnt_type=5:  d2Qdyy
                qnt_type=10: FFQ -> Forcing function

                The following are meant for post-processing:
                -------------------------------------------
                qnt_type=11: FFQ_inv     -> inviscid component of forcing function (divergence of the inviscid flux)
                qnt_type=12: FFQ_vis     -> viscous component of forcing function  (divergence of the viscous flux)
                qnt_type=13: FFQ_SA_prod -> Production term of SA source term
                qnt_type=14: FFQ_SA_dest -> Destruction term of SA source term
                qnt_type=15: FFQ_SA_dist -> Distribution term of SA source term
                qnt_type=16: FFQ_SA_cond -> Conservation term of SA source term

           Argument, qnt_id :           
                           Primitive variable           |  PDE forcing function 
                        or its derivatives              |
                           (qnt_type<=5)                |  (qnt_type>=10)
                      ------------------------------------------------------------------------------------------
                qnt_id=0:  RHO                          |  continuity 
                qnt_id=1:  U                            |  x-momentum
                qnt_id=2:  V                            |  y-momentum
                qnt_id=3:  E                            |  energy
                qnt_id=4:  NTL (nu_tilde)               |  SA
                qnt_id=9:  P
                qnt_id=10: d   (Manufcrd wall distance)
                
      Examples of calls:

           dEdy = MS_qnt_2D(1, 2, 3, 0.1, 0.2); 
           or:
           dEdy = MS_qnt_2D(1, qt_dQdy, qi_E, 0.1, 0.2); //if using suggested enum definitions
           Interpretation: provide me with the first derivative wrt to 'y' of the manufactured field 'E' of MS-1 at the point (0.1,0.2) of the domain.

           Q_SA = MS_qnt_2D(2, 10, 4, 0.5, 0.3);
           or:           
           Q_SA = MS_qnt_2D(2, qt_FFQ, qi_SA, 0.5, 0.3);  //if using suggested enum definitions           
           Interpretation: provide me with the forcing function of MS-2 for the SA equation at the point (0.5,0.3) of the domain.


      Conservative variables: The conservative variables can be computed via primitive variables as RhoE = RHO*E, etc.
                              The derivatives of the conservative variables (to impose the BC for viscous fluxes, etc.) 
                              can be computed by chain rule operations on the derivatives of the primitive variables, ex:
                                                   d(RhoU)dx = d(RHO)dx*U + RHO*d(U)dx

*/


  double qnt;

  switch(MS_no)
  {
  
/********************************************************************************/
/*                              FREE FLOWS                                      */
/********************************************************************************/
/*                                                                              */
/*                              Euler CASES                                     */
/*                                                                              */
/*    MS-1: EULER Verification case subsonic                                    */
/**/case 1:  qnt = MS_1(qnt_type,qnt_id,x,y); break;                          /**/
/*                                                                              */
/*    MS-2: EULER Verification case supersonic on curved domain                 */
/**/case 2:  qnt = MS_2(qnt_type,qnt_id,x,y); break;                          /**/
/*                                                                              */
/*                                                                              */
/*                          Navier-Stokes CASES                                 */
/*                                                                              */
/*    MS-3: NS Verification case                                                */
/**/case 3:  qnt = MS_3(qnt_type,qnt_id,x,y); break;                          /**/
/*                                                                              */    
/*                                                                              */                
/*                            RANS-SA CASES                                     */ 
/*                                                                              */
/*    MS-4: RANS-SA Verification case, original SA model                        */
/**/case 4:  qnt = MS_4(qnt_type,qnt_id,x,y); break;                          /**/
/*                                                                              */
/*    MS-5: RANS-SA Verification case, negative SA model                        */
/**/case 5:  qnt = MS_5(qnt_type,qnt_id,x,y); break;                          /**/
/*                                                                              */
/********************************************************************************/


/********************************************************************************/
/*                           WALL_BOUNDED FLOWS                                 */
/********************************************************************************/
/*                                                                              */
/*                              Euler CASES                                     */
/*                                                                              */
/*    MS-6:  EULER Verification case, slip wall-bounded curved domain           */
/**/case 6:  qnt = MS_6(qnt_type,qnt_id,x,y); break;                          /**/
/*                                                                              */    
/*                                                                              */    
/*                          Navier-Stokes CASES                                 */    
/*                                                                              */    
/*    MS-7:  Navier-Stokes Verification case, no-slip wall-bounded curved domain*/
/**/case 7:  qnt = MS_7(qnt_type,qnt_id,x,y); break;                          /**/
/*                                                                              */
/*                                                                              */
/*                            RANS-SA CASES                                     */ 
/*                                                                              */
/*    MS-8:  Realistic RANS-SA Verification case, Eça et al. (2007)             */
/**/case 8:  qnt = MS_8(qnt_type,qnt_id,x,y); break;                          /**/
/*                                                                              */
/*    MS-9:  Realistic RANS-SA Verification case, Oliver et al. (2012)          */
/**/case 9:  qnt = MS_9(qnt_type,qnt_id,x,y); break;                            */
/*                                                                              */
/********************************************************************************/

    // if nothing is chosen, exit the program.
    default:    printf("No MS is chosen! exit...\n"); exit(1);    
  } 
  
  return(qnt);
  
} 

/********************************************************************************/
/*                              FREE FLOWS                                      */
/********************************************************************************/

//                              Euler CASES                                     //
double MS_1(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution forcing functions and quantities      
      
    // MS-1: EULER Verification case subsonic
      
*/
                        
  const double rho_0   =  1.0;
  const double rho_x   =  0.3;
  const double rho_y   = -0.2;
  const double rho_xy  =  0.3;  
  const double u_0     =  1.0;
  const double u_x     =  0.3;
  const double u_y     =  0.3;
  const double u_xy    =  0.3;
  const double v_0     =  1.0;
  const double v_x     =  0.3;
  const double v_y     =  0.3;
  const double v_xy    =  0.3;  
  const double p_0     =  18.0;
  const double p_x     =  5.0; 
  const double p_y     =  5.0; 
  const double p_xy    =  0.5;
 
  const double a_rhox  =  1.0;
  const double a_rhoy  =  1.0;
  const double a_rhoxy =  1.0;  
  const double a_ux    =  3.0;
  const double a_uy    =  1.0;
  const double a_uxy   =  1.0;
  const double a_vx    =  1.0;
  const double a_vy    =  1.0;
  const double a_vxy   =  1.0;  
  const double a_px    =  2.0;
  const double a_py    =  1.0;
  const double a_pxy   =  1.0; 
  
  const double L       =  1.0;

// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.
  const double Gamma   =  1.4;
  
  double RHO,U,V,E,P;
  double dRHOdx,dRHOdy;
  double dUdx,dUdy;
  double dVdx,dVdy;  
  double dEdx,dEdy; 
  double dPdx,dPdy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E;

  double qnt;
  int MS_no;
  
  MS_no = 1;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = rho_0 + rho_x*sin(piMS*a_rhox*x/L) + rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + rho_y*cos(piMS*a_rhoy*y/L);
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = piMS*a_rhox*rho_x*cos(piMS*a_rhox*x/L)/L - piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L)/L; 
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dRHOdy = -piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*y/L)*cos(piMS*a_rhoxy*x/L)/L - piMS*a_rhoy*rho_y*sin(piMS*a_rhoy*y/L)/L;
          qnt    = dRHOdy;
          break;      

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break;                    
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = u_0 + u_x*sin(piMS*a_ux*x/L) + u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L) + u_y*cos(piMS*a_uy*y/L); 
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = piMS*a_ux*u_x*cos(piMS*a_ux*x/L)/L - piMS*a_uxy*u_xy*sin(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L)/L;           
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dUdy = -piMS*a_uxy*u_xy*sin(piMS*a_uxy*y/L)*cos(piMS*a_uxy*x/L)/L - piMS*a_uy*u_y*sin(piMS*a_uy*y/L)/L; 
          qnt  = dUdy;
          break;                     

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);     

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    
          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U;
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_U_VIS = 0.;
          
          qnt    = FFQ_U_VIS;
          break;              
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V = v_0 + v_x*cos(piMS*a_vx*x/L) + v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L) + v_y*sin(piMS*a_vy*y/L);
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -piMS*a_vx*v_x*sin(piMS*a_vx*x/L)/L - piMS*a_vxy*v_xy*sin(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L)/L;
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = -piMS*a_vxy*v_xy*sin(piMS*a_vxy*y/L)*cos(piMS*a_vxy*x/L)/L + piMS*a_vy*v_y*cos(piMS*a_vy*y/L)/L;
          qnt  = dVdy;
          break;                   
         
      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);   
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V; 
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_V_VIS = 0.; 
          qnt    = FFQ_V_VIS;
          break;               
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
      
        
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO);
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO);           
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1)*pow(RHO, 2)) + dPdy/((Gamma - 1)*RHO);          
          qnt = dEdy;
          break;      
   
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;          
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_E_VIS = 0;
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;      
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          P = p_0 + p_x*cos(piMS*a_px*x/L) + p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L) + p_y*sin(piMS*a_py*y/L); 
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
          dPdx = -piMS*a_px*p_x*sin(piMS*a_px*x/L)/L - piMS*a_pxy*p_xy*sin(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L)/L; 
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = -piMS*a_pxy*p_xy*sin(piMS*a_pxy*y/L)*cos(piMS*a_pxy*x/L)/L + piMS*a_py*p_y*cos(piMS*a_py*y/L)/L; 
          qnt  = dPdy;
          break;          

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;            
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}

double MS_2(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution forcing functions and quantities      
      
    // MS-2: EULER Verification case supersonic on curved domain
      
*/
  const double rho_0   =  2.7; 
  const double rho_x   =  0.9; 
  const double rho_y   = -0.9; 
  const double rho_xy  =  1.0;                         
  const double u_0     =  2.0; 
  const double u_x     =  0.7; 
  const double u_y     =  0.7; 
  const double u_xy    =  0.4; 
  const double v_0     =  2.0; 
  const double v_x     =  0.4; 
  const double v_y     =  0.4; 
  const double v_xy    =  0.4; 
  const double p_0     =  2.0; 
  const double p_x     =  1.0; 
  const double p_y     =  1.0; 
  const double p_xy    =  0.5; 

  const double a_rhox  =  1.5; 
  const double a_rhoy  =  1.5; 
  const double a_rhoxy =  1.5; 
  const double a_ux    =  1.0; 
  const double a_uy    =  1.0; 
  const double a_uxy   =  1.0; 
  const double a_vx    =  1.0; 
  const double a_vy    =  1.0; 
  const double a_vxy   =  1.0; 
  const double a_px    =  1.0; 
  const double a_py    =  1.0; 
  const double a_pxy   =  1.5; 
  
  const double L       =  1.0;

// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.
  const double Gamma   =  1.4;  
    
  double RHO,U,V,E,P;
  double dRHOdx,dRHOdy;
  double dUdx,dUdy;
  double dVdx,dVdy;  
  double dEdx,dEdy; 
  double dPdx,dPdy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E;

  double qnt;
  int MS_no;
  
  MS_no = 2;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = rho_0 + rho_x*sin(piMS*a_rhox*x/L) + rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + rho_y*cos(piMS*a_rhoy*y/L);
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = piMS*a_rhox*rho_x*cos(piMS*a_rhox*x/L)/L - piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L)/L; 
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dRHOdy = -piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*y/L)*cos(piMS*a_rhoxy*x/L)/L - piMS*a_rhoy*rho_y*sin(piMS*a_rhoy*y/L)/L; 
          qnt    = dRHOdy;
          break;      

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break;                    
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = u_0 + u_x*sin(piMS*a_ux*x/L) + u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L) + u_y*cos(piMS*a_uy*y/L); 
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = piMS*a_ux*u_x*cos(piMS*a_ux*x/L)/L - piMS*a_uxy*u_xy*sin(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L)/L; 
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dUdy = -piMS*a_uxy*u_xy*sin(piMS*a_uxy*y/L)*cos(piMS*a_uxy*x/L)/L - piMS*a_uy*u_y*sin(piMS*a_uy*y/L)/L; 
          qnt  = dUdy;
          break;                     

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);     

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    
          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U;          
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_U_VIS = 0.;
          
          qnt    = FFQ_U_VIS;
          break;              
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V = v_0 + v_x*cos(piMS*a_vx*x/L) + v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L) + v_y*sin(piMS*a_vy*y/L); 
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -piMS*a_vx*v_x*sin(piMS*a_vx*x/L)/L - piMS*a_vxy*v_xy*sin(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L)/L; 
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = -piMS*a_vxy*v_xy*sin(piMS*a_vxy*y/L)*cos(piMS*a_vxy*x/L)/L + piMS*a_vy*v_y*cos(piMS*a_vy*y/L)/L;
          qnt  = dVdy;
          break;                   
         
      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);   
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V; 
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_V_VIS = 0.; 
          qnt    = FFQ_V_VIS;
          break;               
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
      
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO); 
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO); 
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdy/((Gamma - 1.0)*RHO); 
          qnt = dEdy;
          break;      
   
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_E_VIS = 0;
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;
      
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          P = p_0 + p_x*cos(piMS*a_px*x/L) + p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L) + p_y*sin(piMS*a_py*y/L);
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
          dPdx = -piMS*a_px*p_x*sin(piMS*a_px*x/L)/L - piMS*a_pxy*p_xy*sin(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L)/L; 
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = -piMS*a_pxy*p_xy*sin(piMS*a_pxy*y/L)*cos(piMS*a_pxy*x/L)/L + piMS*a_py*p_y*cos(piMS*a_py*y/L)/L; 
          qnt  = dPdy;
          break;          

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;            
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}


//                          Navier-Stokes CASES                                 //   

double MS_3(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution forcing functions and quantities      
      
    // MS-3 NS Verification case
      
*/

  const double rho_0   =  1.0;
  const double rho_x   =  0.1;
  const double rho_y   = -0.2;
  const double rho_xy  =  0.1;  
  const double u_0     =  2.0; 
  const double u_x     =  0.3;
  const double u_y     =  0.3; 
  const double u_xy    =  0.3;
  const double v_0     =  2.0;
  const double v_x     =  0.3;
  const double v_y     =  0.3; 
  const double v_xy    =  0.3;  
  const double p_0     =  10.0;
  const double p_x     =  1.0; 
  const double p_y     =  1.0; 
  const double p_xy    =  0.5;
  
  const double a_rhox  =  1.0;
  const double a_rhoy  =  1.0;
  const double a_rhoxy =  1.0;  
  const double a_ux    =  3.0;
  const double a_uy    =  1.0;
  const double a_uxy   =  1.0;
  const double a_vx    =  1.0;
  const double a_vy    =  1.0;
  const double a_vxy   =  1.0;
  const double a_px    =  2.0;
  const double a_py    =  1.0;
  const double a_pxy   =  1.0;    

  const double L       =  1.0;

// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.
  const double Gamma   =  1.4;
  const double Pr      =  0.7;  
  const double mu      =  1e-1;        
  
  
  double RHO,U,V,E,P;
  double dRHOdx,dRHOdy,d2RHOdxx,d2RHOdyy;
  double dUdx,dUdy,d2Udxx,d2Udxy,d2Udyy;
  double dVdx,dVdy,d2Vdxx,d2Vdxy,d2Vdyy;   
  double dEdx,dEdy,d2Edxx,d2Edyy; // double d2Edxy; 
  double dPdx,dPdy,d2Pdxx,d2Pdyy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E;

  double qnt;
  int MS_no;
  
  MS_no = 3;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = rho_0 + rho_x*sin(piMS*a_rhox*x/L) + rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + rho_y*cos(piMS*a_rhoy*y/L);
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = piMS*a_rhox*rho_x*cos(piMS*a_rhox*x/L)/L - piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L)/L; 
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dRHOdy = -piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*y/L)*cos(piMS*a_rhoxy*x/L)/L - piMS*a_rhoy*rho_y*sin(piMS*a_rhoy*y/L)/L; 
          qnt    = dRHOdy;
          break; 

        // qnt_type=3: d2Qdxx
        case 3:    
          d2RHOdxx = -pow(piMS, 2)*(pow(a_rhox, 2)*rho_x*sin(piMS*a_rhox*x/L) + pow(a_rhoxy, 2)*rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L))/pow(L, 2);
          qnt = d2RHOdxx;
          break;

        // qnt_type=4: d2Qdxy
//         case 4:    
//           d2RHOdxy = ;
//           qnt = d2RHOdxy;
//           break;
          
        // qnt_type=5: d2Qdyy
        case 5:    
          d2RHOdyy = -pow(piMS, 2)*(pow(a_rhoxy, 2)*rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + pow(a_rhoy, 2)*rho_y*cos(piMS*a_rhoy*y/L))/pow(L, 2); 
          qnt = d2RHOdyy;
          break;           

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break;                    
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = u_0 + u_x*sin(piMS*a_ux*x/L) + u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L) + u_y*cos(piMS*a_uy*y/L); 
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = piMS*a_ux*u_x*cos(piMS*a_ux*x/L)/L - piMS*a_uxy*u_xy*sin(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L)/L; 
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dUdy = -piMS*a_uxy*u_xy*sin(piMS*a_uxy*y/L)*cos(piMS*a_uxy*x/L)/L - piMS*a_uy*u_y*sin(piMS*a_uy*y/L)/L; 
          qnt  = dUdy;
          break;   
          
        // qnt_type=3: d2Qdxx
        case 3:              
          d2Udxx = -pow(piMS, 2)*(pow(a_ux, 2)*u_x*sin(piMS*a_ux*x/L) + pow(a_uxy, 2)*u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L))/pow(L, 2); 
          qnt    = d2Udxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:    
          d2Udxy = pow(piMS, 2)*pow(a_uxy, 2)*u_xy*sin(piMS*a_uxy*x/L)*sin(piMS*a_uxy*y/L)/pow(L, 2);
          qnt    = d2Udxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:    
          d2Udyy = -pow(piMS, 2)*(pow(a_uxy, 2)*u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L) + pow(a_uy, 2)*u_y*cos(piMS*a_uy*y/L))/pow(L, 2); 
          qnt    = d2Udyy;
          break;           

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);     

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    
          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U;          
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);          
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);          

          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);            

          FFQ_U_VIS = -1.0L/3.0L*mu*(4*d2Udxx + 3*d2Udyy + d2Vdxy);
          qnt    = FFQ_U_VIS;
          break;                        
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V = v_0 + v_x*cos(piMS*a_vx*x/L) + v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L) + v_y*sin(piMS*a_vy*y/L); 
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -piMS*a_vx*v_x*sin(piMS*a_vx*x/L)/L - piMS*a_vxy*v_xy*sin(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L)/L; 
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = -piMS*a_vxy*v_xy*sin(piMS*a_vxy*y/L)*cos(piMS*a_vxy*x/L)/L + piMS*a_vy*v_y*cos(piMS*a_vy*y/L)/L;
          qnt  = dVdy;
          break;
          
        // qnt_type=3: d2Qdxx
        case 3:    
          d2Vdxx = -pow(piMS, 2)*(pow(a_vx, 2)*v_x*cos(piMS*a_vx*x/L) + pow(a_vxy, 2)*v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L))/pow(L, 2);
          qnt    = d2Vdxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:    
          d2Vdxy = pow(piMS, 2)*pow(a_vxy, 2)*v_xy*sin(piMS*a_vxy*x/L)*sin(piMS*a_vxy*y/L)/pow(L, 2);
          qnt    = d2Vdxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2Vdyy = -pow(piMS, 2)*(pow(a_vxy, 2)*v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L) + pow(a_vy, 2)*v_y*sin(piMS*a_vy*y/L))/pow(L, 2);
          qnt    = d2Vdyy;
          break;            
         
      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);   
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2.*dVdy*RHO*V; 
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);          
  
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);                          
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);                               

          FFQ_V_VIS = -1.0L/3.0L*mu*(d2Udxy + 3*d2Vdxx + 4*d2Vdyy); 
          qnt    = FFQ_V_VIS;
          break;               
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
      
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO); 
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO); 
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdy/((Gamma - 1.0)*RHO);  
          qnt = dEdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:  
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Rho,x,y);
          d2RHOdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Rho,x,y);
          
          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          d2Udxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y); 
          d2Vdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);            
                   
          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);
          dPdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_Prsr,x,y);
          d2Pdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Prsr,x,y); 
                    
          d2Edxx = 1.0*U*d2Udxx + 1.0*V*d2Vdxx + 1.0*pow(dUdx, 2.0) + 1.0*pow(dVdx, 2.0) - P*d2RHOdxx/((Gamma - 1.0)*pow(RHO, 2.0)) + 2.0*P*pow(dRHOdx, 2.0)/((Gamma - 1.0)*pow(RHO, 3.0)) + d2Pdxx/((Gamma - 1.0)*RHO) - 2.0*dPdx*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0));
          qnt = d2Edxx;
          break;

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Rho,x,y);
          d2RHOdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Rho,x,y);          

          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);      
          dPdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_Prsr,x,y);  
          d2Pdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Prsr,x,y);                       
                   
          d2Edyy = 1.0*U*d2Udyy + 1.0*V*d2Vdyy + 1.0*pow(dUdy, 2.0) + 1.0*pow(dVdy, 2.0) - P*d2RHOdyy/((Gamma - 1.0)*pow(RHO, 2.0)) + 2.0*P*pow(dRHOdy, 2.0)/((Gamma - 1.0)*pow(RHO, 3.0)) + d2Pdyy/((Gamma - 1.0)*RHO) - 2.0*dPdy*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)); 
          qnt    = d2Edyy;
          break;                
   
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          U      = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          d2Edxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_E,x,y);
          d2Edyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_E,x,y);                
                   
          FFQ_E_VIS = (1.0L/3.0L)*mu*(Pr*(-2.*dUdx*(2.*dUdx - dVdy) - 3.*dUdy*(dUdy + dVdx) - 3.*dVdx*(dUdy + dVdx) + 2.*dVdy*(dUdx - 2.*dVdy) - 2.*(2.*d2Udxx - d2Vdxy)*U - 3.*(d2Udxy + d2Vdxx)*V + 2.*(d2Udxy - 2*d2Vdyy)*V - 3.*(d2Udyy + d2Vdxy)*U) + 3.*Gamma*(-d2Edxx + d2Udxx*U + d2Vdxx*V + pow(dUdx, 2) + pow(dVdx, 2)) + 3.*Gamma*(-d2Edyy + d2Udyy*U + d2Vdyy*V + pow(dUdy, 2) + pow(dVdy, 2)))/Pr; 
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;
      
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          P = p_0 + p_x*cos(piMS*a_px*x/L) + p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L) + p_y*sin(piMS*a_py*y/L);
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
          dPdx = -piMS*a_px*p_x*sin(piMS*a_px*x/L)/L - piMS*a_pxy*p_xy*sin(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L)/L; 
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = -piMS*a_pxy*p_xy*sin(piMS*a_pxy*y/L)*cos(piMS*a_pxy*x/L)/L + piMS*a_py*p_y*cos(piMS*a_py*y/L)/L; 
          qnt  = dPdy;
          break; 
          
        // qnt_type=3: d2Qdxx
        case 3:    
          d2Pdxx = -pow(piMS, 2)*(pow(a_px, 2)*p_x*cos(piMS*a_px*x/L) + pow(a_pxy, 2)*p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L))/pow(L, 2); 
          qnt = d2Pdxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:                        
          d2Pdyy = -pow(piMS, 2)*(pow(a_pxy, 2)*p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L) + pow(a_py, 2)*p_y*sin(piMS*a_py*y/L))/pow(L, 2);
          qnt    = d2Pdyy;
          break;          

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;            
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}


//                            RANS-SA CASES                                    // 
double MS_4(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution forcing functions and quantities      
      
    // MS-4: RANS-SA Verification case, original SA model
      
*/  
  const double rho_0    =  1.0;
  const double rho_x    =  0.1;
  const double rho_y    = -0.2;
  const double rho_xy   =  0.1;  
  const double u_0      =  2.0; 
  const double u_x      =  0.3;
  const double u_y      =  0.3; 
  const double u_xy     =  0.3;
  const double v_0      =  2.0;
  const double v_x      =  0.3;
  const double v_y      =  0.3; 
  const double v_xy     =  0.3;  
  const double p_0      =  10.0;
  const double p_x      =  1.0; 
  const double p_y      =  1.0; 
  const double p_xy     =  0.5;
  const double nu_sa_0  =  6.0e-1;
  const double nu_sa_x  = -0.3e-1;
  const double nu_sa_y  = -0.2e-1;
  const double nu_sa_xy =  0.2e-1;  
  
  
  const double a_rhox   =  1.0;
  const double a_rhoy   =  1.0;
  const double a_rhoxy  =  1.0;  
  const double a_ux     =  3.0;
  const double a_uy     =  1.0;
  const double a_uxy    =  1.0;
  const double a_vx     =  1.0;
  const double a_vy     =  1.0;
  const double a_vxy    =  1.0;
  const double a_px     =  2.0;
  const double a_py     =  1.0;
  const double a_pxy    =  1.0;   
  const double a_nusax  =  2.0;
  const double a_nusay  =  1.0;  
  const double a_nusaxy =  3.0;     

  const double L        =  1.0;

// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.                        
  const double Gamma    =  1.4;
  const double Pr       =  0.7;    
  const double mu       =  1e-3;
                      
// RANS-SA variables
  const double sigma    = 2./3.;
  const double Prt      = 0.9  ;  
  const double c_b1     = 0.1355;
  const double c_b2     = 0.622;
  const double c_t3     = 1.2;
  const double c_t4     = 0.5;       
  const double c_v1     = 7.1;
  const double kappa    = 0.41;     
  const double c_w1     = c_b1/(kappa*kappa)+(1.0+c_b2)/sigma;
  const double c_w2     = 0.3;
  const double c_w3     = 2.0;  
  const double r_lim    = 10.0;
                        
  const double d0       = 1.e0;    // fixed wall distance
  
  double RHO,U,V,E,P,NTL;
  double d;
  double dRHOdx,dRHOdy,d2RHOdxx,d2RHOdyy;  
  double dUdx,dUdy,d2Udxx,d2Udxy,d2Udyy;
  double dVdx,dVdy,d2Vdxx,d2Vdxy,d2Vdyy;  
  double dEdx,dEdy,d2Edxx,d2Edyy; // double d2Edxy; 
  double dNTLdx,dNTLdy,d2NTLdxx,d2NTLdyy; // double d2NTLdxy;
  double dPdx,dPdy,d2Pdxx,d2Pdyy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV,FFQ_NTL_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS,FFQ_NTL_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E,FFQ_NTL;
  double FFQ_NTL_SRC_PRODUCT,FFQ_NTL_SRC_DESTRUCT,FFQ_NTL_SRC_DISTRIB,FFQ_NTL_SRC_CONSERV;   
  
  double qnt;
  int MS_no;
  
  MS_no = 4;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = rho_0 + rho_x*sin(piMS*a_rhox*x/L) + rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + rho_y*cos(piMS*a_rhoy*y/L); 
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = piMS*a_rhox*rho_x*cos(piMS*a_rhox*x/L)/L - piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L)/L; 
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dRHOdy = -piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*y/L)*cos(piMS*a_rhoxy*x/L)/L - piMS*a_rhoy*rho_y*sin(piMS*a_rhoy*y/L)/L;
          qnt    = dRHOdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:
          d2RHOdxx = -pow(piMS, 2)*(pow(a_rhox, 2)*rho_x*sin(piMS*a_rhox*x/L) + pow(a_rhoxy, 2)*rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L))/pow(L, 2);
          qnt    = d2RHOdxx;
          break;      

        // qnt_type=4: d2Qdxy
//        case 4:    
//          qnt    = d2RHOdxy;
//          break;
          
        // qnt_type=5: d2Qdyy
        case 5:
          d2RHOdyy = -pow(piMS, 2)*(pow(a_rhoxy, 2)*rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + pow(a_rhoy, 2)*rho_y*cos(piMS*a_rhoy*y/L))/pow(L, 2); 
          qnt    = d2RHOdyy;
          break;           

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);          
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break; 
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = u_0 + u_x*sin(piMS*a_ux*x/L) + u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L) + u_y*cos(piMS*a_uy*y/L); 
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = piMS*a_ux*u_x*cos(piMS*a_ux*x/L)/L - piMS*a_uxy*u_xy*sin(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L)/L; 
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dUdy = -piMS*a_uxy*u_xy*sin(piMS*a_uxy*y/L)*cos(piMS*a_uxy*x/L)/L - piMS*a_uy*u_y*sin(piMS*a_uy*y/L)/L;
          qnt  = dUdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2Udxx = -pow(piMS, 2)*(pow(a_ux, 2)*u_x*sin(piMS*a_ux*x/L) + pow(a_uxy, 2)*u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L))/pow(L, 2); 
          qnt    = d2Udxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:    
          d2Udxy = pow(piMS, 2)*pow(a_uxy, 2)*u_xy*sin(piMS*a_uxy*x/L)*sin(piMS*a_uxy*y/L)/pow(L, 2); 
          qnt    = d2Udxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:    
          d2Udyy = -pow(piMS, 2)*(pow(a_uxy, 2)*u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L) + pow(a_uy, 2)*u_y*cos(piMS*a_uy*y/L))/pow(L, 2);
          qnt    = d2Udyy;
          break;                

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);      

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                      
                    
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U; 
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);          
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);          

          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);  
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);            
                 
          FFQ_U_VIS = ((2.0L/3.0L)*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(-(2*d2Udxx - d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) - (2*dUdx - dVdy)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(-(d2Udyy + d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) - (dUdy + dVdx)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(2*(2*dUdx - dVdy)*(dNTLdx*RHO + dRHOdx*NTL) + 3*(dUdy + dVdx)*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2))/pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2);           
          qnt    = FFQ_U_VIS;
          break;                          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V = v_0 + v_x*cos(piMS*a_vx*x/L) + v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L) + v_y*sin(piMS*a_vy*y/L);
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -piMS*a_vx*v_x*sin(piMS*a_vx*x/L)/L - piMS*a_vxy*v_xy*sin(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L)/L;
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = -piMS*a_vxy*v_xy*sin(piMS*a_vxy*y/L)*cos(piMS*a_vxy*x/L)/L + piMS*a_vy*v_y*cos(piMS*a_vy*y/L)/L; 
          qnt  = dVdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2Vdxx = -pow(piMS, 2)*(pow(a_vx, 2)*v_x*cos(piMS*a_vx*x/L) + pow(a_vxy, 2)*v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L))/pow(L, 2); 
          qnt    = d2Vdxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:              
          d2Vdxy = pow(piMS, 2)*pow(a_vxy, 2)*v_xy*sin(piMS*a_vxy*x/L)*sin(piMS*a_vxy*y/L)/pow(L, 2);
          qnt    = d2Vdxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2Vdyy = -pow(piMS, 2)*(pow(a_vxy, 2)*v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L) + pow(a_vy, 2)*v_y*sin(piMS*a_vy*y/L))/pow(L, 2);
          qnt    = d2Vdyy;
          break;                

      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);             
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);             
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V; 
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);          

          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);   
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);                          
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y); 
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);
          
          FFQ_V_VIS = (-(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udxy + d2Vdxx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdy + dVdx)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (2.0L/3.0L)*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udxy - 2*d2Vdyy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdx - 2*dVdy)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(-2*(dUdx - 2*dVdy)*(dNTLdy*RHO + dRHOdy*NTL) + 3*(dUdy + dVdx)*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2))/pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2); 
          
          qnt    = FFQ_V_VIS;
          break;          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO); 
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO); 
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdy/((Gamma - 1.0)*RHO); 
          qnt = dEdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:  
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Rho,x,y);
          d2RHOdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Rho,x,y);
          
          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          d2Udxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y); 
          d2Vdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);            
                   
          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);
          dPdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_Prsr,x,y);
          d2Pdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Prsr,x,y); 
                    
          d2Edxx = 1.0*U*d2Udxx + 1.0*V*d2Vdxx + 1.0*pow(dUdx, 2.0) + 1.0*pow(dVdx, 2.0) - P*d2RHOdxx/((Gamma - 1.0)*pow(RHO, 2.0)) + 2.0*P*pow(dRHOdx, 2.0)/((Gamma - 1.0)*pow(RHO, 3.0)) + d2Pdxx/((Gamma - 1.0)*RHO) - 2.0*dPdx*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0));
          qnt = d2Edxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Rho,x,y);
          d2RHOdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Rho,x,y);          

          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);      
          dPdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_Prsr,x,y);  
          d2Pdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Prsr,x,y);                       
                   
          d2Edyy = 1.0*U*d2Udyy + 1.0*V*d2Vdyy + 1.0*pow(dUdy, 2.0) + 1.0*pow(dVdy, 2.0) - P*d2RHOdyy/((Gamma - 1.0)*pow(RHO, 2.0)) + 2.0*P*pow(dRHOdy, 2.0)/((Gamma - 1.0)*pow(RHO, 3.0)) + d2Pdyy/((Gamma - 1.0)*RHO) - 2.0*dPdy*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)); 
          qnt    = d2Edyy;
          break;        
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          E      = MS_qnt_2D(MS_no,qt_Q,     qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_E,x,y);
          d2Edxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_E,x,y);
          d2Edyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_E,x,y);
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y); 
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                      
          
          
          FFQ_E_VIS = ((2.0L/3.0L)*Pr*Prt*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(-dUdx*(2*dUdx - dVdy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + dVdy*(dUdx - 2*dVdy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) - (2*d2Udxx - d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*U + (d2Udxy - 2*d2Vdyy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*V + (dUdx - 2*dVdy)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)*V - (2*dUdx - dVdy)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)*U) - Pr*Prt*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(dUdy*(dUdy + dVdx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + dVdx*(dUdy + dVdx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (d2Udxy + d2Vdxx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*V + (d2Udyy + d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*U + (dUdy + dVdx)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)*V + (dUdy + dVdx)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)*U) + Pr*Prt*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(-2*(dUdx - 2*dVdy)*(dNTLdy*RHO + dRHOdy*NTL)*V + 2*(2*dUdx - dVdy)*(dNTLdx*RHO + dRHOdx*NTL)*U + 3*(dUdy + dVdx)*(dNTLdx*RHO + dRHOdx*NTL)*V + 3*(dUdy + dVdx)*(dNTLdy*RHO + dRHOdy*NTL)*U)*pow(NTL, 2)*pow(RHO, 2) + Gamma*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*(-d2Edxx + 1.0*d2Udxx*U + 1.0*d2Vdxx*V + 1.0*pow(dUdx, 2) + 1.0*pow(dVdx, 2)) + (Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*(-d2Edyy + 1.0*d2Udyy*U + 1.0*d2Vdyy*V + 1.0*pow(dUdy, 2) + 1.0*pow(dVdy, 2)) + (-dEdx + 1.0*dUdx*U + 1.0*dVdx*V)*(4*Pr*dNTLdx*NTL*pow(RHO, 2) + 4*Pr*dRHOdx*pow(NTL, 2)*RHO + 3*Prt*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2) + (-dEdy + 1.0*dUdy*U + 1.0*dVdy*V)*(4*Pr*dNTLdy*NTL*pow(RHO, 2) + 4*Pr*dRHOdy*pow(NTL, 2)*RHO + 3*Prt*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) - 3*Gamma*((dNTLdx*RHO + dRHOdx*NTL)*(-dEdx + 1.0*dUdx*U + 1.0*dVdx*V) + (dNTLdy*RHO + dRHOdy*NTL)*(-dEdy + 1.0*dUdy*U + 1.0*dVdy*V))*(Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*pow(NTL, 2)*pow(RHO, 2))/(Pr*Prt*pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2));
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;
      
    // qnt_id=4: NTL
    case 4:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          NTL = nu_sa_0 + nu_sa_x*cos(piMS*a_nusax*x/L) + nu_sa_xy*cos(piMS*a_nusaxy*x/L)*cos(piMS*a_nusaxy*y/L) + nu_sa_y*cos(piMS*a_nusay*y/L); 
          qnt = NTL;
          break;
          
        // qnt_type=1: dQdx
        case 1:                        
          dNTLdx = -piMS*a_nusax*nu_sa_x*sin(piMS*a_nusax*x/L)/L - piMS*a_nusaxy*nu_sa_xy*sin(piMS*a_nusaxy*x/L)*cos(piMS*a_nusaxy*y/L)/L;
          qnt  = dNTLdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dNTLdy = -piMS*a_nusaxy*nu_sa_xy*sin(piMS*a_nusaxy*y/L)*cos(piMS*a_nusaxy*x/L)/L - piMS*a_nusay*nu_sa_y*sin(piMS*a_nusay*y/L)/L; 
          qnt  = dNTLdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2NTLdxx = -pow(piMS, 2)*(pow(a_nusax, 2)*nu_sa_x*cos(piMS*a_nusax*x/L) + pow(a_nusaxy, 2)*nu_sa_xy*cos(piMS*a_nusaxy*x/L)*cos(piMS*a_nusaxy*y/L))/pow(L, 2);
          qnt    = d2NTLdxx;
          break;      

        // qnt_type=4: d2Qdxy
//         case 4:    
//           d2NTLdxy = 4*sig_v*y*(pow(sig_v, 2)*pow(y, 2)/pow(x, 2) - 1)*exp(-pow(sig_v, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 3)); 
//           qnt    = d2NTLdxy;
//           break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2NTLdyy = -pow(piMS, 2)*(pow(a_nusaxy, 2)*nu_sa_xy*cos(piMS*a_nusaxy*x/L)*cos(piMS*a_nusaxy*y/L) + pow(a_nusay, 2)*nu_sa_y*cos(piMS*a_nusay*y/L))/pow(L, 2); 
          qnt    = d2NTLdyy;
          break;                

     // qnt_type=10: forcing function
        case 10:

          FFQ_NTL_INV          = MS_qnt_2D(MS_no,qt_FFQ_inv,    qi_SA,x,y);
          FFQ_NTL_VIS          = MS_qnt_2D(MS_no,qt_FFQ_vis,    qi_SA,x,y);
          FFQ_NTL_SRC_PRODUCT  = MS_qnt_2D(MS_no,qt_FFQ_SA_prod,qi_SA,x,y);
          FFQ_NTL_SRC_DESTRUCT = MS_qnt_2D(MS_no,qt_FFQ_SA_dest,qi_SA,x,y);
          FFQ_NTL_SRC_DISTRIB  = MS_qnt_2D(MS_no,qt_FFQ_SA_dist,qi_SA,x,y);
          FFQ_NTL_SRC_CONSERV  = MS_qnt_2D(MS_no,qt_FFQ_SA_cons,qi_SA,x,y);
          
          FFQ_NTL = FFQ_NTL_INV + FFQ_NTL_VIS + FFQ_NTL_SRC_PRODUCT  + FFQ_NTL_SRC_DESTRUCT + FFQ_NTL_SRC_DISTRIB + FFQ_NTL_SRC_CONSERV; 
          qnt    = FFQ_NTL;
          break;
       
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U        = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          
          V        = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                       

          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                

          FFQ_NTL_INV = dNTLdx*RHO*U + dNTLdy*RHO*V + dRHOdx*NTL*U + dRHOdy*NTL*V + dUdx*NTL*RHO + dVdy*NTL*RHO; 
          qnt    = FFQ_NTL_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
        
          NTL      = MS_qnt_2D(MS_no,qt_Q,     qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Ntl,x,y);                       
          d2NTLdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Ntl,x,y);
          d2NTLdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Ntl,x,y);
                    
          FFQ_NTL_VIS = -(d2NTLdxx*(mu + NTL*RHO) + d2NTLdyy*(mu + NTL*RHO) + dNTLdx*(dNTLdx*RHO + dRHOdx*NTL) + dNTLdy*(dNTLdy*RHO + dRHOdy*NTL))/sigma;          
          qnt = FFQ_NTL_VIS;
          break;
          
        // qnt_type=13: Production term in the SA source term
        case 13:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);          
          
          FFQ_NTL_SRC_PRODUCT = -c_b1*(-c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2)) + 1)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*NTL*RHO;         
          qnt = FFQ_NTL_SRC_PRODUCT;
          break;       
          
        // qnt_type=14: Destruction term in the SA source term
        case 14:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);           
          
          FFQ_NTL_SRC_DESTRUCT = (-c_b1*c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2))/pow(kappa, 2) + c_w1*pow((pow(c_w3, 6) + 1)/(pow(c_w3, 6) + pow(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6)), 1.0L/6.0L)*(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))))*pow(NTL, 2)*RHO/pow(d, 2);           
          qnt = FFQ_NTL_SRC_DESTRUCT;
          break;  

        // qnt_type=15: Distribtuion term in the SA source term 
        case 15:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                 
          
          FFQ_NTL_SRC_DISTRIB = -c_b2*(pow(dNTLdx, 2) + pow(dNTLdy, 2))*RHO/sigma;           
          qnt = FFQ_NTL_SRC_DISTRIB;
          break;            
          
        // qnt_type=13: Conservation term in the SA source term 
        case 16:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
        
          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                       
                    
          FFQ_NTL_SRC_CONSERV = (dNTLdx*dRHOdx + dNTLdy*dRHOdy)*(mu/RHO + NTL)/sigma;
          qnt = FFQ_NTL_SRC_CONSERV;
          break;               
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;       
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
       case 0:              
          P = p_0 + p_x*cos(piMS*a_px*x/L) + p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L) + p_y*sin(piMS*a_py*y/L); 
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
	  dPdx = -piMS*a_px*p_x*sin(piMS*a_px*x/L)/L - piMS*a_pxy*p_xy*sin(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L)/L;
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = -piMS*a_pxy*p_xy*sin(piMS*a_pxy*y/L)*cos(piMS*a_pxy*x/L)/L + piMS*a_py*p_y*cos(piMS*a_py*y/L)/L;
          qnt  = dPdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:    
          d2Pdxx = -pow(piMS, 2)*(pow(a_px, 2)*p_x*cos(piMS*a_px*x/L) + pow(a_pxy, 2)*p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L))/pow(L, 2);
          qnt = d2Pdxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:                        
          d2Pdyy = -pow(piMS, 2)*(pow(a_pxy, 2)*p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L) + pow(a_py, 2)*p_y*sin(piMS*a_py*y/L))/pow(L, 2);
          qnt    = d2Pdyy;
          break;               

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;  
      
    // qnt_id=10: d
    case 10:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          d   = y+d0;
          qnt = d;
          break;
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;        
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}

double MS_5(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution forcing functions and quantities      
      
    // MS-5: RANS-SA Verification case, modified SA model
      
*/  
  const double rho_0    =  1.0;
  const double rho_x    =  0.1;
  const double rho_y    = -0.2;
  const double rho_xy   =  0.1;  
  const double u_0      =  2.0; 
  const double u_x      =  0.3;
  const double u_y      =  0.3; 
  const double u_xy     =  0.3;
  const double v_0      =  2.0;
  const double v_x      =  0.3;
  const double v_y      =  0.3; 
  const double v_xy     =  0.3;  
  const double p_0      =  10.0;
  const double p_x      =  1.0; 
  const double p_y      =  1.0; 
  const double p_xy     =  0.5;
  const double nu_sa_0  = -6.0e-0;
  const double nu_sa_x  = -0.3e-0;
  const double nu_sa_y  = -0.2e-0;
  const double nu_sa_xy =  0.2e-0;  
  
  
  const double a_rhox   =  1.0;
  const double a_rhoy   =  1.0;
  const double a_rhoxy  =  1.0;  
  const double a_ux     =  3.0;
  const double a_uy     =  1.0;
  const double a_uxy    =  1.0;
  const double a_vx     =  1.0;
  const double a_vy     =  1.0;
  const double a_vxy    =  1.0;
  const double a_px     =  2.0;
  const double a_py     =  1.0;
  const double a_pxy    =  1.0;  
  const double a_nusax  =  2.0;  
  const double a_nusay  =  1.0;
  const double a_nusaxy =  3.0;      
                        
  const double L        =  1.0;
        
// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.                        
  const double Gamma    =  1.4;
  const double Pr       =  0.7; 
  const double mu       =  1e-1;
                        
// RANS-SA variables    
  const double sigma    = 2./3.;
  const double c_b1     = 0.1355;
  const double c_b2     = 0.622;
  const double c_t3     = 1.2;
  const double kappa    = 0.41;     
  const double c_w1     = c_b1/(kappa*kappa)+(1.0+c_b2)/sigma;
  const double c_n1     = 16.;  
                        
  const double d0       = 1.e0;    // fixed wall distance
  
  double RHO,U,V,E,P,NTL;
  double d;
  double dRHOdx,dRHOdy,d2RHOdxx,d2RHOdyy;  
  double dUdx,dUdy,d2Udxx,d2Udxy,d2Udyy;
  double dVdx,dVdy,d2Vdxx,d2Vdxy,d2Vdyy;  
  double dEdx,dEdy,d2Edxx,d2Edyy; // double d2Edxy; 
  double dNTLdx,dNTLdy,d2NTLdxx,d2NTLdyy; // double d2NTLdxy; 
  double dPdx,dPdy,d2Pdxx,d2Pdyy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV,FFQ_NTL_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS,FFQ_NTL_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E,FFQ_NTL;  
  double FFQ_NTL_SRC_PRODUCT,FFQ_NTL_SRC_DESTRUCT,FFQ_NTL_SRC_DISTRIB,FFQ_NTL_SRC_CONSERV; 
  
  double qnt;
  int MS_no;
  
  MS_no = 5;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = rho_0 + rho_x*sin(piMS*a_rhox*x/L) + rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + rho_y*cos(piMS*a_rhoy*y/L); 
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = piMS*a_rhox*rho_x*cos(piMS*a_rhox*x/L)/L - piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L)/L; 
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dRHOdy = -piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*y/L)*cos(piMS*a_rhoxy*x/L)/L - piMS*a_rhoy*rho_y*sin(piMS*a_rhoy*y/L)/L;
          qnt    = dRHOdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:
          d2RHOdxx = -pow(piMS, 2)*(pow(a_rhox, 2)*rho_x*sin(piMS*a_rhox*x/L) + pow(a_rhoxy, 2)*rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L))/pow(L, 2);
          qnt    = d2RHOdxx;
          break;      

        // qnt_type=4: d2Qdxy
//        case 4:    
//          qnt    = d2RHOdxy;
//          break;
          
        // qnt_type=5: d2Qdyy
        case 5:
          d2RHOdyy = -pow(piMS, 2)*(pow(a_rhoxy, 2)*rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + pow(a_rhoy, 2)*rho_y*cos(piMS*a_rhoy*y/L))/pow(L, 2); 
          qnt    = d2RHOdyy;
          break;           

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);          
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break; 
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = u_0 + u_x*sin(piMS*a_ux*x/L) + u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L) + u_y*cos(piMS*a_uy*y/L); 
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = piMS*a_ux*u_x*cos(piMS*a_ux*x/L)/L - piMS*a_uxy*u_xy*sin(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L)/L; 
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dUdy = -piMS*a_uxy*u_xy*sin(piMS*a_uxy*y/L)*cos(piMS*a_uxy*x/L)/L - piMS*a_uy*u_y*sin(piMS*a_uy*y/L)/L;
          qnt  = dUdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2Udxx = -pow(piMS, 2)*(pow(a_ux, 2)*u_x*sin(piMS*a_ux*x/L) + pow(a_uxy, 2)*u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L))/pow(L, 2); 
          qnt    = d2Udxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:    
          d2Udxy = pow(piMS, 2)*pow(a_uxy, 2)*u_xy*sin(piMS*a_uxy*x/L)*sin(piMS*a_uxy*y/L)/pow(L, 2); 
          qnt    = d2Udxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:    
          d2Udyy = -pow(piMS, 2)*(pow(a_uxy, 2)*u_xy*cos(piMS*a_uxy*x/L)*cos(piMS*a_uxy*y/L) + pow(a_uy, 2)*u_y*cos(piMS*a_uy*y/L))/pow(L, 2);
          qnt    = d2Udyy;
          break;                

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);      

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                      
                    
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U; 
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);          
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);          
  
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);         
                 
          FFQ_U_VIS = -1.0L/3.0L*mu*(4*d2Udxx + 3*d2Udyy + d2Vdxy);
          qnt    = FFQ_U_VIS;
          break;                          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V = v_0 + v_x*cos(piMS*a_vx*x/L) + v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L) + v_y*sin(piMS*a_vy*y/L);
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -piMS*a_vx*v_x*sin(piMS*a_vx*x/L)/L - piMS*a_vxy*v_xy*sin(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L)/L;
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = -piMS*a_vxy*v_xy*sin(piMS*a_vxy*y/L)*cos(piMS*a_vxy*x/L)/L + piMS*a_vy*v_y*cos(piMS*a_vy*y/L)/L; 
          qnt  = dVdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2Vdxx = -pow(piMS, 2)*(pow(a_vx, 2)*v_x*cos(piMS*a_vx*x/L) + pow(a_vxy, 2)*v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L))/pow(L, 2); 
          qnt    = d2Vdxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:              
          d2Vdxy = pow(piMS, 2)*pow(a_vxy, 2)*v_xy*sin(piMS*a_vxy*x/L)*sin(piMS*a_vxy*y/L)/pow(L, 2);
          qnt    = d2Vdxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2Vdyy = -pow(piMS, 2)*(pow(a_vxy, 2)*v_xy*cos(piMS*a_vxy*x/L)*cos(piMS*a_vxy*y/L) + pow(a_vy, 2)*v_y*sin(piMS*a_vy*y/L))/pow(L, 2);
          qnt    = d2Vdyy;
          break;                

      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);             
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);             
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V;
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);          

          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);                          
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);
                  
          FFQ_V_VIS = -1.0L/3.0L*mu*(d2Udxy + 3*d2Vdxx + 4*d2Vdyy);
          qnt    = FFQ_V_VIS;
          break;          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO); 
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO); 
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdy/((Gamma - 1.0)*RHO); 
          qnt = dEdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:  
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Rho,x,y);
          d2RHOdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Rho,x,y);
          
          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          d2Udxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y); 
          d2Vdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);            
                   
          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);
          dPdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_Prsr,x,y);
          d2Pdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Prsr,x,y); 
                    
          d2Edxx = 1.0*U*d2Udxx + 1.0*V*d2Vdxx + 1.0*pow(dUdx, 2) + 1.0*pow(dVdx, 2) - P*d2RHOdxx/((Gamma - 1)*pow(RHO, 2)) + 2*P*pow(dRHOdx, 2)/((Gamma - 1)*pow(RHO, 3)) + d2Pdxx/((Gamma - 1)*RHO) - 2*dPdx*dRHOdx/((Gamma - 1)*pow(RHO, 2)); 
          qnt = d2Edxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Rho,x,y);
          d2RHOdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Rho,x,y);          

          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);      
          dPdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_Prsr,x,y);  
          d2Pdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Prsr,x,y);                       
                   
          d2Edyy = 1.0*U*d2Udyy + 1.0*V*d2Vdyy + 1.0*pow(dUdy, 2) + 1.0*pow(dVdy, 2) - P*d2RHOdyy/((Gamma - 1)*pow(RHO, 2)) + 2*P*pow(dRHOdy, 2)/((Gamma - 1)*pow(RHO, 3)) + d2Pdyy/((Gamma - 1)*RHO) - 2*dPdy*dRHOdy/((Gamma - 1)*pow(RHO, 2));
          qnt    = d2Edyy;
          break;        
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;          
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          E      = MS_qnt_2D(MS_no,qt_Q,     qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_E,x,y);
          d2Edxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_E,x,y);
          d2Edyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_E,x,y);
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y); 
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                      
          
          
          FFQ_E_VIS = (1.0L/3.0L)*mu*(Pr*(-2*dUdx*(2*dUdx - dVdy) - 3*dUdy*(dUdy + dVdx) - 3*dVdx*(dUdy + dVdx) + 2*dVdy*(dUdx - 2*dVdy) - 2*(2*d2Udxx - d2Vdxy)*U - 3*(d2Udxy + d2Vdxx)*V + 2*(d2Udxy - 2*d2Vdyy)*V - 3*(d2Udyy + d2Vdxy)*U) + 3*Gamma*(-d2Edxx + d2Udxx*U + d2Vdxx*V + pow(dUdx, 2) + pow(dVdx, 2)) + 3*Gamma*(-d2Edyy + d2Udyy*U + d2Vdyy*V + pow(dUdy, 2) + pow(dVdy, 2)))/Pr; 
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;
      
    // qnt_id=4: NTL
    case 4:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          NTL = nu_sa_0 + nu_sa_x*cos(piMS*a_nusax*x/L) + nu_sa_xy*cos(piMS*a_nusaxy*x/L)*cos(piMS*a_nusaxy*y/L) + nu_sa_y*cos(piMS*a_nusay*y/L); 
          qnt = NTL;
          break;
          
        // qnt_type=1: dQdx
        case 1:                        
          dNTLdx = -piMS*a_nusax*nu_sa_x*sin(piMS*a_nusax*x/L)/L - piMS*a_nusaxy*nu_sa_xy*sin(piMS*a_nusaxy*x/L)*cos(piMS*a_nusaxy*y/L)/L;  
          qnt  = dNTLdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dNTLdy = -piMS*a_nusaxy*nu_sa_xy*sin(piMS*a_nusaxy*y/L)*cos(piMS*a_nusaxy*x/L)/L - piMS*a_nusay*nu_sa_y*sin(piMS*a_nusay*y/L)/L;
          qnt  = dNTLdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2NTLdxx = -pow(piMS, 2)*(pow(a_nusax, 2)*nu_sa_x*cos(piMS*a_nusax*x/L) + pow(a_nusaxy, 2)*nu_sa_xy*cos(piMS*a_nusaxy*x/L)*cos(piMS*a_nusaxy*y/L))/pow(L, 2);
          qnt    = d2NTLdxx;
          break;      

        // qnt_type=4: d2Qdxy
//         case 4:    
//           d2NTLdxy = 4*sig_v*y*(pow(sig_v, 2)*pow(y, 2)/pow(x, 2) - 1)*exp(-pow(sig_v, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 3)); 
//           qnt    = d2NTLdxy;
//           break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2NTLdyy = -pow(piMS, 2)*(pow(a_nusaxy, 2)*nu_sa_xy*cos(piMS*a_nusaxy*x/L)*cos(piMS*a_nusaxy*y/L) + pow(a_nusay, 2)*nu_sa_y*cos(piMS*a_nusay*y/L))/pow(L, 2); 
          qnt    = d2NTLdyy;
          break;                

    // qnt_type=10: forcing function
        case 10:

          FFQ_NTL_INV          = MS_qnt_2D(MS_no,qt_FFQ_inv,    qi_SA,x,y);
          FFQ_NTL_VIS          = MS_qnt_2D(MS_no,qt_FFQ_vis,    qi_SA,x,y);
          FFQ_NTL_SRC_PRODUCT  = MS_qnt_2D(MS_no,qt_FFQ_SA_prod,qi_SA,x,y);
          FFQ_NTL_SRC_DESTRUCT = MS_qnt_2D(MS_no,qt_FFQ_SA_dest,qi_SA,x,y);
          FFQ_NTL_SRC_DISTRIB  = MS_qnt_2D(MS_no,qt_FFQ_SA_dist,qi_SA,x,y);
          FFQ_NTL_SRC_CONSERV  = MS_qnt_2D(MS_no,qt_FFQ_SA_cons,qi_SA,x,y);
          
          FFQ_NTL = FFQ_NTL_INV + FFQ_NTL_VIS + FFQ_NTL_SRC_PRODUCT  + FFQ_NTL_SRC_DESTRUCT + FFQ_NTL_SRC_DISTRIB + FFQ_NTL_SRC_CONSERV; 
          qnt    = FFQ_NTL;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U        = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          
          V        = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                       

          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                

          FFQ_NTL_INV = dNTLdx*RHO*U + dNTLdy*RHO*V + dRHOdx*NTL*U + dRHOdy*NTL*V + dUdx*NTL*RHO + dVdy*NTL*RHO; 
          qnt    = FFQ_NTL_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,     qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Ntl,x,y);                       
          d2NTLdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Ntl,x,y);
          d2NTLdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Ntl,x,y);
          
          FFQ_NTL_VIS = -(dNTLdx*((c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3))*(dNTLdx*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*RHO + dRHOdx*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*NTL + 3*(dNTLdx*RHO + dRHOdx*NTL)*pow(NTL, 3)*pow(RHO, 3)) + 3*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(dNTLdx*RHO + dRHOdx*NTL)*pow(NTL, 3)*pow(RHO, 3)) + dNTLdy*((c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3))*(dNTLdy*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*RHO + dRHOdy*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*NTL + 3*(dNTLdy*RHO + dRHOdy*NTL)*pow(NTL, 3)*pow(RHO, 3)) + 3*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(dNTLdy*RHO + dRHOdy*NTL)*pow(NTL, 3)*pow(RHO, 3)) + (d2NTLdxx + d2NTLdyy)*(c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3))*(mu*(c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3)) + (c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*NTL*RHO))/(sigma*pow(c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3), 2)); 
          qnt = FFQ_NTL_VIS;
          break;
          
        // qnt_type=13: Production term in the SA source term
        case 13:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);          
          
          FFQ_NTL_SRC_PRODUCT = -c_b1*(-c_t3 + 1)*NTL*RHO*fabs(dUdy - dVdx);  
          qnt = FFQ_NTL_SRC_PRODUCT;
          break;       
          
        // qnt_type=14: Destruction term in the SA source term
        case 14:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);            
          
          FFQ_NTL_SRC_DESTRUCT = -c_w1*pow(NTL, 2)*RHO/pow(d, 2); 
          qnt = FFQ_NTL_SRC_DESTRUCT;
          break;  

        // qnt_type=15: Distribtuion term in the SA source term 
        case 15:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                       
          
          FFQ_NTL_SRC_DISTRIB = -c_b2*(pow(dNTLdx, 2) + pow(dNTLdy, 2))*RHO/sigma;
          qnt = FFQ_NTL_SRC_DISTRIB;
          break;            
          
        // qnt_type=13: Conservation term in the SA source term 
        case 16:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
          
          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                       
          
          
          FFQ_NTL_SRC_CONSERV = (dNTLdx*dRHOdx + dNTLdy*dRHOdy)*(mu/RHO + (c_n1 + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))*NTL/(c_n1 - pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3)))/sigma; 
          qnt = FFQ_NTL_SRC_CONSERV;
          break;               
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;       
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
       case 0:              
          P = p_0 + p_x*cos(piMS*a_px*x/L) + p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L) + p_y*sin(piMS*a_py*y/L); 
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
	  dPdx = -piMS*a_px*p_x*sin(piMS*a_px*x/L)/L - piMS*a_pxy*p_xy*sin(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L)/L;
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = -piMS*a_pxy*p_xy*sin(piMS*a_pxy*y/L)*cos(piMS*a_pxy*x/L)/L + piMS*a_py*p_y*cos(piMS*a_py*y/L)/L;
          qnt  = dPdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:    
          d2Pdxx = -pow(piMS, 2)*(pow(a_px, 2)*p_x*cos(piMS*a_px*x/L) + pow(a_pxy, 2)*p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L))/pow(L, 2); 
          qnt = d2Pdxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:                        
          d2Pdyy = -pow(piMS, 2)*(pow(a_pxy, 2)*p_xy*cos(piMS*a_pxy*x/L)*cos(piMS*a_pxy*y/L) + pow(a_py, 2)*p_y*sin(piMS*a_py*y/L))/pow(L, 2);
          qnt    = d2Pdyy;
          break;               

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;  
      
    // qnt_id=10: d
    case 10:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          d   = y+d0;
          qnt = d;
          break;
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;        
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}

/********************************************************************************/
/*                           WALL_BOUNDED FLOWS                                 */
/********************************************************************************/

//                              Euler CASES                                     //
double MS_6(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution and forcing functions      
      
    // MS-6: EULER Verification case, slip wall-bounded curved domain
      
*/
  

// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.  
  const double Gamma   =  1.4;


  const double rho_0   =  1.0;  
  const double p_0     =  1.0;

  // Constants used in the spatial deformation functions
  double const A_def = 0.05;
  double const N_def = 2.00;
  double const L_def = 1.00;
  
  double const u_wall = 1.0; // To ensure non-zero wall velocity in the inviscid mode or when userdef is used.

  
  double RHO,U,V,E,P;
  double dRHOdx,dRHOdy;
  double dUdx,dUdy;
  double dVdx,dVdy;  
  double dEdx,dEdy; 
  double dPdx,dPdy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E;

  double qnt;
  int MS_no;
  
  MS_no = 6;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = rho_0 + pow(-A_def*sin(piMS*N_def*x/L_def) + y, 2); 
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = -2*piMS*A_def*N_def*(-A_def*sin(piMS*N_def*x/L_def) + y)*cos(piMS*N_def*x/L_def)/L_def;
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:              
          dRHOdy = -2*A_def*sin(piMS*N_def*x/L_def) + 2*y; 
          qnt    = dRHOdy;
          break;      

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO;
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break;                    
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = -A_def*sin(piMS*N_def*x/L_def) + u_wall + y; 
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = -piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def;
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dUdy = 1;
          qnt  = dUdy;
          break;                     

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);     

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    
          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U;
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_U_VIS = 0.;
          
          qnt    = FFQ_U_VIS;
          break;              
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V = piMS*A_def*N_def*(-A_def*sin(piMS*N_def*x/L_def) + u_wall + y)*cos(piMS*N_def*x/L_def)/L_def; 
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -pow(piMS, 2)*pow(A_def, 2)*pow(N_def, 2)*pow(cos(piMS*N_def*x/L_def), 2)/pow(L_def, 2) - pow(piMS, 2)*A_def*pow(N_def, 2)*(-A_def*sin(piMS*N_def*x/L_def) + u_wall + y)*sin(piMS*N_def*x/L_def)/pow(L_def, 2); 
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def;
          qnt  = dVdy;
          break;                   
         
      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);   
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V;
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_V_VIS = 0.; 
          qnt    = FFQ_V_VIS;
          break;               
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
      
        
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO);
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO);           
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1)*pow(RHO, 2)) + dPdy/((Gamma - 1)*RHO);          
          qnt = dEdy;
          break;      
   
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          FFQ_E_VIS = 0;
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;      
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          P = p_0 + pow(-A_def*sin(piMS*N_def*x/L_def) + y, 2);
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
          dPdx = -2*piMS*A_def*N_def*(-A_def*sin(piMS*N_def*x/L_def) + y)*cos(piMS*N_def*x/L_def)/L_def;
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = -2*A_def*sin(piMS*N_def*x/L_def) + 2*y;
          qnt  = dPdy;
          break;          

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;            
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}


//                          Navier-Stokes CASES                                 //   

double MS_7(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution and forcing functions      
      
    // MS-7: Navier-Stokes Verification case, no-slip wall-bounded curved domain
      
*/

// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.
  const double Gamma   =  1.4;
  const double Pr      =  0.7;  
  const double mu      =  1e-5; 

  const double rho_0   =  1.2;  
  const double rho_c   =  3.0;   
  const double p_0     =  2.0;  // To stabilizer the pressure (which is otherwise 0 in outer BL regions) on coarse grids.
  
  const double sig_u   =  2.0;   

  // Constants used in the spatial deformation functions
  double const A_def = 0.05;
  double const N_def = 2.00;
  double const L_def = 1.00;


  
  double RHO,U,V,E,P;
  double dRHOdx,dRHOdy,d2RHOdxx,d2RHOdyy;
  double dUdx,dUdy,d2Udxx,d2Udxy,d2Udyy;
  double dVdx,dVdy,d2Vdxx,d2Vdxy,d2Vdyy;   
  double dEdx,dEdy,d2Edxx,d2Edyy; // double d2Edxy; 
  double dPdx,dPdy,d2Pdxx,d2Pdyy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E;

  double qnt;
  int MS_no;
  
  MS_no = 7;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = rho_0 + pow(-A_def*sin(piMS*N_def*x/L_def) + y, 2)/pow(rho_c, 2); 
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = -2*piMS*A_def*N_def*(-A_def*sin(piMS*N_def*x/L_def) + y)*cos(piMS*N_def*x/L_def)/(L_def*pow(rho_c, 2));
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dRHOdy = (-2*A_def*sin(piMS*N_def*x/L_def) + 2*y)/pow(rho_c, 2);
          qnt    = dRHOdy;
          break; 

        // qnt_type=3: d2Qdxx
        case 3:    
          d2RHOdxx = 2*pow(piMS, 2)*A_def*pow(N_def, 2)*(A_def*pow(cos(piMS*N_def*x/L_def), 2) - (A_def*sin(piMS*N_def*x/L_def) - y)*sin(piMS*N_def*x/L_def))/(pow(L_def, 2)*pow(rho_c, 2)); 
          qnt = d2RHOdxx;
          break;

        // qnt_type=4: d2Qdxy
//         case 4:    
//           d2RHOdxy = ;
//           qnt = d2RHOdxy;
//           break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2RHOdyy = 2/pow(rho_c, 2);
          qnt = d2RHOdyy;
          break;           

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break;                    
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = -erf(1.0*sig_u*(A_def*sin(piMS*N_def*x/L_def) - y)/sqrt(x));
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = -2*(1.0*piMS*A_def*N_def*sig_u*cos(piMS*N_def*x/L_def)/(L_def*sqrt(x)) - 0.5*sig_u*(A_def*sin(piMS*N_def*x/L_def) - y)/pow(x, 3.0L/2.0L))*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)/sqrt(piMS);
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dUdy = 2.0*sig_u*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)/(sqrt(piMS)*sqrt(x));
          qnt  = dUdy;
          break;   
          
        // qnt_type=3: d2Qdxx
        case 3:              
          d2Udxx = 2*sig_u*(1.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/(L_def*x) + 1.0*pow(piMS, 2)*A_def*pow(N_def, 2)*sin(piMS*N_def*x/L_def)/pow(L_def, 2) + pow(sig_u, 2)*(A_def*sin(piMS*N_def*x/L_def) - y)*(1.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def - (0.5*A_def*sin(piMS*N_def*x/L_def) - 0.5*y)/x)*(2.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def - (1.0*A_def*sin(piMS*N_def*x/L_def) - 1.0*y)/x)/x - (0.75*A_def*sin(piMS*N_def*x/L_def) - 0.75*y)/pow(x, 2))*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)/(sqrt(piMS)*sqrt(x)); 
          qnt    = d2Udxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:    
          d2Udxy = -sig_u*(4.0*pow(sig_u, 2)*(A_def*sin(piMS*N_def*x/L_def) - y)*(1.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def - (0.5*A_def*sin(piMS*N_def*x/L_def) - 0.5*y)/x) + 1.0)*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)/(sqrt(piMS)*pow(x, 3.0L/2.0L)); 
          qnt    = d2Udxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:    
          d2Udyy = 4.0*pow(sig_u, 3)*(A_def*sin(piMS*N_def*x/L_def) - y)*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)/(sqrt(piMS)*pow(x, 3.0L/2.0L)); 
          qnt    = d2Udyy;
          break;           

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);     

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    
          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U;          
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);          
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);          

          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);            

          FFQ_U_VIS = -1.0L/3.0L*mu*(4*d2Udxx + 3*d2Udyy + d2Vdxy);
          qnt    = FFQ_U_VIS;
          break;                        
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V = -piMS*A_def*N_def*cos(piMS*N_def*x/L_def)*erf(1.0*sig_u*(A_def*sin(piMS*N_def*x/L_def) - y)/sqrt(x))/L_def; 
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -2*sqrt(piMS)*A_def*N_def*(1.0*piMS*A_def*N_def*sig_u*cos(piMS*N_def*x/L_def)/(L_def*sqrt(x)) - 0.5*sig_u*(A_def*sin(piMS*N_def*x/L_def) - y)/pow(x, 3.0L/2.0L))*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)*cos(piMS*N_def*x/L_def)/L_def + pow(piMS, 2)*A_def*pow(N_def, 2)*sin(piMS*N_def*x/L_def)*erf(1.0*sig_u*(A_def*sin(piMS*N_def*x/L_def) - y)/sqrt(x))/pow(L_def, 2); 
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = 2.0*sqrt(piMS)*A_def*N_def*sig_u*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)*cos(piMS*N_def*x/L_def)/(L_def*sqrt(x)); 
          qnt  = dVdy;
          break;
          
        // qnt_type=3: d2Qdxx
        case 3:    
          d2Vdxx = A_def*N_def*(2*sqrt(piMS)*pow(sig_u, 3)*(A_def*sin(piMS*N_def*x/L_def) - y)*(1.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def - (0.5*A_def*sin(piMS*N_def*x/L_def) - 0.5*y)/x)*(2.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def - (1.0*A_def*sin(piMS*N_def*x/L_def) - 1.0*y)/x)*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)*cos(piMS*N_def*x/L_def)/pow(x, 3.0L/2.0L) + 2*sqrt(piMS)*sig_u*(1.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/(L_def*x) + 1.0*pow(piMS, 2)*A_def*pow(N_def, 2)*sin(piMS*N_def*x/L_def)/pow(L_def, 2) - (0.75*A_def*sin(piMS*N_def*x/L_def) - 0.75*y)/pow(x, 2))*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)*cos(piMS*N_def*x/L_def)/sqrt(x) + 4*pow(piMS, 3.0L/2.0L)*N_def*sig_u*(1.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def - (0.5*A_def*sin(piMS*N_def*x/L_def) - 0.5*y)/x)*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)*sin(piMS*N_def*x/L_def)/(L_def*sqrt(x)) + pow(piMS, 3)*pow(N_def, 2)*cos(piMS*N_def*x/L_def)*erf(1.0*sig_u*(A_def*sin(piMS*N_def*x/L_def) - y)/sqrt(x))/pow(L_def, 2))/L_def; 
          qnt    = d2Vdxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:    
          d2Vdxy = -sqrt(piMS)*A_def*N_def*sig_u*(4.0*pow(sig_u, 2)*(A_def*sin(piMS*N_def*x/L_def) - y)*(1.0*piMS*A_def*N_def*cos(piMS*N_def*x/L_def)/L_def - (0.5*A_def*sin(piMS*N_def*x/L_def) - 0.5*y)/x)*cos(piMS*N_def*x/L_def)/x + 1.0*cos(piMS*N_def*x/L_def)/x + 2.0*piMS*N_def*sin(piMS*N_def*x/L_def)/L_def)*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)/(L_def*sqrt(x));
          qnt    = d2Vdxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2Vdyy = 4.0*sqrt(piMS)*A_def*N_def*pow(sig_u, 3)*(A_def*sin(piMS*N_def*x/L_def) - y)*exp(-1.0*pow(sig_u, 2)*pow(A_def*sin(piMS*N_def*x/L_def) - y, 2)/x)*cos(piMS*N_def*x/L_def)/(L_def*pow(x, 3.0L/2.0L));
          qnt    = d2Vdyy;
          break;            
         
      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);   
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2.*dVdy*RHO*V; 
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);          
  
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);                          
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);                               

          FFQ_V_VIS = -1.0L/3.0L*mu*(d2Udxy + 3*d2Vdxx + 4*d2Vdyy); 
          qnt    = FFQ_V_VIS;
          break;               
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
      
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO); 
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO); 
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdy/((Gamma - 1.0)*RHO);  
          qnt = dEdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:  
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Rho,x,y);
          d2RHOdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Rho,x,y);
          
          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          d2Udxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y); 
          d2Vdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);            
                   
          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);
          dPdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_Prsr,x,y);
          d2Pdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Prsr,x,y); 
                    
          d2Edxx = 1.0*U*d2Udxx + 1.0*V*d2Vdxx + 1.0*pow(dUdx, 2.0) + 1.0*pow(dVdx, 2.0) - P*d2RHOdxx/((Gamma - 1.0)*pow(RHO, 2.0)) + 2.0*P*pow(dRHOdx, 2.0)/((Gamma - 1.0)*pow(RHO, 3.0)) + d2Pdxx/((Gamma - 1.0)*RHO) - 2.0*dPdx*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0));
          qnt = d2Edxx;
          break;

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Rho,x,y);
          d2RHOdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Rho,x,y);          

          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);      
          dPdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_Prsr,x,y);  
          d2Pdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Prsr,x,y);                       
                   
          d2Edyy = 1.0*U*d2Udyy + 1.0*V*d2Vdyy + 1.0*pow(dUdy, 2.0) + 1.0*pow(dVdy, 2.0) - P*d2RHOdyy/((Gamma - 1.0)*pow(RHO, 2.0)) + 2.0*P*pow(dRHOdy, 2.0)/((Gamma - 1.0)*pow(RHO, 3.0)) + d2Pdyy/((Gamma - 1.0)*RHO) - 2.0*dPdy*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)); 
          qnt    = d2Edyy;
          break;                
   
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          U      = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          d2Edxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_E,x,y);
          d2Edyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_E,x,y);                
                   
          FFQ_E_VIS = (1.0L/3.0L)*mu*(Pr*(-2.*dUdx*(2.*dUdx - dVdy) - 3.*dUdy*(dUdy + dVdx) - 3.*dVdx*(dUdy + dVdx) + 2.*dVdy*(dUdx - 2.*dVdy) - 2.*(2.*d2Udxx - d2Vdxy)*U - 3.*(d2Udxy + d2Vdxx)*V + 2.*(d2Udxy - 2*d2Vdyy)*V - 3.*(d2Udyy + d2Vdxy)*U) + 3.*Gamma*(-d2Edxx + d2Udxx*U + d2Vdxx*V + pow(dUdx, 2) + pow(dVdx, 2)) + 3.*Gamma*(-d2Edyy + d2Udyy*U + d2Vdyy*V + pow(dUdy, 2) + pow(dVdy, 2)))/Pr; 
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;
      
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          P = p_0 + pow(-A_def*sin(piMS*N_def*x/L_def) + y, 2); 
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
          dPdx = -2*piMS*A_def*N_def*(-A_def*sin(piMS*N_def*x/L_def) + y)*cos(piMS*N_def*x/L_def)/L_def;
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = -2*A_def*sin(piMS*N_def*x/L_def) + 2*y; 
          qnt  = dPdy;
          break; 
          
        // qnt_type=3: d2Qdxx
        case 3:    
          d2Pdxx = 2*pow(piMS, 2)*A_def*pow(N_def, 2)*(A_def*pow(cos(piMS*N_def*x/L_def), 2) - (A_def*sin(piMS*N_def*x/L_def) - y)*sin(piMS*N_def*x/L_def))/pow(L_def, 2); 
          qnt = d2Pdxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:                        
          d2Pdyy = 2;
          qnt    = d2Pdyy;
          break;          

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;            
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}


//                            RANS-SA CASES                                    // 
double MS_8(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution and forcing functions      
      
    // MS-8: Realistic RANS-SA Verification case, Eça et al. (2007)     
      
*/  

                        
  const double rho_0   =  1.1;
  const double rho_x   =  0.2;
  const double rho_y   =  0.5;
  const double rho_xy  =  0.2;  
  const double a_rhox  =  1.0;
  const double a_rhoy  =  2.0;
  const double a_rhoxy =  3.0; 
  const double L       =  1.0;  

// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.
  const double Gamma   =  1.4;
  const double Pr      =  0.7;
  const double mu      =  1.e-6;  
  
  const double sig_u   =  4.0;  
  
  const double p_0     =  0.1;  // To stabilizer the pressure on coarse grids.  
  
  const double sig_nu  =  2.5*sig_u;   
  const double NTL_max =  1.e-3;   
    
  // RANS-SA variables
  const double sigma   = 2./3.;
  const double Prt     = 0.9  ;  
  const double c_b1    = 0.1355;
  const double c_b2    = 0.622;
  const double c_t3    = 1.2;
  const double c_t4    = 0.5;       
  const double c_v1    = 7.1;
  const double c_v2    = 0.7;
  const double c_v3    = 0.9;
  const double kappa   = 0.41;     
  const double c_w1    = c_b1/(kappa*kappa)+(1.0+c_b2)/sigma;
  const double c_w2    = 0.3;
  const double c_w3    = 2.0;  
  const double r_lim   = 10.0;
  
  double S_MS,Sbar_MS,nu_MS,chi_MS,fv1_MS,fv2_MS;    
  double RHO,U,V,E,P,NTL;
  double d;
  double dRHOdx,dRHOdy,d2RHOdxx,d2RHOdyy;  
  double dUdx,dUdy,d2Udxx,d2Udxy,d2Udyy;
  double dVdx,dVdy,d2Vdxx,d2Vdxy,d2Vdyy;  
  double dEdx,dEdy,d2Edxx,d2Edyy; // double d2Edxy;
  double dNTLdx,dNTLdy,d2NTLdxx,d2NTLdyy; // double d2NTLdxy; 
  double dPdx,dPdy,d2Pdxx,d2Pdyy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV,FFQ_NTL_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS,FFQ_NTL_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E,FFQ_NTL;
  double FFQ_NTL_SRC_PRODUCT,FFQ_NTL_SRC_DESTRUCT,FFQ_NTL_SRC_DISTRIB,FFQ_NTL_SRC_CONSERV;   
  double tau_w;  
  
  double qnt;
  int MS_no;
  
  MS_no = 8;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = rho_0 + rho_x*sin(piMS*a_rhox*x/L) + rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + rho_y*cos(piMS*a_rhoy*y/L);  
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = piMS*a_rhox*rho_x*cos(piMS*a_rhox*x/L)/L - piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L)/L; 
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:              
          dRHOdy = -piMS*a_rhoxy*rho_xy*sin(piMS*a_rhoxy*y/L)*cos(piMS*a_rhoxy*x/L)/L - piMS*a_rhoy*rho_y*sin(piMS*a_rhoy*y/L)/L;
          qnt    = dRHOdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:
          d2RHOdxx = -pow(piMS, 2)*(pow(a_rhox, 2)*rho_x*sin(piMS*a_rhox*x/L) + pow(a_rhoxy, 2)*rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L))/pow(L, 2);
          qnt    = d2RHOdxx;
          break;      

        // qnt_type=4: d2Qdxy
//        case 4:    
//          qnt    = d2RHOdxy;
//          break;
          
        // qnt_type=5: d2Qdyy
        case 5:          
          d2RHOdyy = -pow(piMS, 2)*(pow(a_rhoxy, 2)*rho_xy*cos(piMS*a_rhoxy*x/L)*cos(piMS*a_rhoxy*y/L) + pow(a_rhoy, 2)*rho_y*cos(piMS*a_rhoy*y/L))/pow(L, 2); 
          qnt    = d2RHOdyy;
          break;           

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);          
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break; 
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = erf(sig_u*y/x); 
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = -2*sig_u*y*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 2));
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:                    
          dUdy = 2*sig_u*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*x); 
          qnt  = dUdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2Udxx = 4*sig_u*y*(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2) + 1)*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 3));
          qnt    = d2Udxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:    
          d2Udxy = 2*sig_u*(2*pow(sig_u, 2)*pow(y, 2)/pow(x, 2) - 1)*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 2));
          qnt    = d2Udxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:    
          d2Udyy = -4*pow(sig_u, 3)*y*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 3)); 
          qnt    = d2Udyy;
          break;                

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);      

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                      
                    
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U; 
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);          
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);          

          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);  
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);            
                 
          FFQ_U_VIS = (-2.0L/3.0L*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((2*d2Udxx - d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (2*dUdx - dVdy)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)) - (pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udyy + d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdy + dVdx)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(2*(2*dUdx - dVdy)*(dNTLdx*RHO + dRHOdx*NTL) + 3*(dUdy + dVdx)*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2))/pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2);           
          qnt    = FFQ_U_VIS;
          break;                          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V = (1 - exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2)))/(sqrt(piMS)*sig_u); 
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -2*sig_u*pow(y, 2)*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 3)); 
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = 2*sig_u*y*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 2)); 
          qnt  = dVdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2Vdxx = 2*sig_u*pow(y, 2)*(-2*pow(sig_u, 2)*pow(y, 2)/pow(x, 2) + 3)*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 4));
          qnt    = d2Vdxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:              
          d2Vdxy = 4*sig_u*y*(pow(sig_u, 2)*pow(y, 2)/pow(x, 2) - 1)*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 3));
          qnt    = d2Vdxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2Vdyy = 2*sig_u*(-2*pow(sig_u, 2)*pow(y, 2)/pow(x, 2) + 1)*exp(-pow(sig_u, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 2));
          qnt    = d2Vdyy;
          break;                

      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);             
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);             
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V; 
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);          

          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);   
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);                          
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y); 
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);
          
          FFQ_V_VIS = (-(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udxy + d2Vdxx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdy + dVdx)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (2.0L/3.0L)*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udxy - 2*d2Vdyy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdx - 2*dVdy)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(-2*(dUdx - 2*dVdy)*(dNTLdy*RHO + dRHOdy*NTL) + 3*(dUdy + dVdx)*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2))/pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2); 
          qnt    = FFQ_V_VIS;
          break;          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO); 
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO); 
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdy/((Gamma - 1.0)*RHO); 
          qnt = dEdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:  
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Rho,x,y);
          d2RHOdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Rho,x,y);
          
          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          d2Udxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y); 
          d2Vdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);            
                   
          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);
          dPdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_Prsr,x,y);
          d2Pdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Prsr,x,y); 
                    
          d2Edxx = 1.0*U*d2Udxx + 1.0*V*d2Vdxx + 1.0*pow(dUdx, 2) + 1.0*pow(dVdx, 2) - P*d2RHOdxx/((Gamma - 1)*pow(RHO, 2)) + 2*P*pow(dRHOdx, 2)/((Gamma - 1)*pow(RHO, 3)) + d2Pdxx/((Gamma - 1)*RHO) - 2*dPdx*dRHOdx/((Gamma - 1)*pow(RHO, 2));
          qnt = d2Edxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Rho,x,y);
          d2RHOdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Rho,x,y);          

          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);      
          dPdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_Prsr,x,y);  
          d2Pdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Prsr,x,y);                       
                   
          d2Edyy = 1.0*U*d2Udyy + 1.0*V*d2Vdyy + 1.0*pow(dUdy, 2) + 1.0*pow(dVdy, 2) - P*d2RHOdyy/((Gamma - 1)*pow(RHO, 2)) + 2*P*pow(dRHOdy, 2)/((Gamma - 1)*pow(RHO, 3)) + d2Pdyy/((Gamma - 1)*RHO) - 2*dPdy*dRHOdy/((Gamma - 1)*pow(RHO, 2)); 
          qnt    = d2Edyy;
          break;        
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          E      = MS_qnt_2D(MS_no,qt_Q,     qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_E,x,y);
          d2Edxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_E,x,y);
          d2Edyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_E,x,y);
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y); 
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                
          
          FFQ_E_VIS = ((2.0L/3.0L)*Pr*Prt*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(-dUdx*(2*dUdx - dVdy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + dVdy*(dUdx - 2*dVdy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) - (2*d2Udxx - d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*U + (d2Udxy - 2*d2Vdyy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*V + (dUdx - 2*dVdy)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)*V - (2*dUdx - dVdy)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)*U) - Pr*Prt*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(dUdy*(dUdy + dVdx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + dVdx*(dUdy + dVdx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (d2Udxy + d2Vdxx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*V + (d2Udyy + d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*U + (dUdy + dVdx)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)*V + (dUdy + dVdx)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)*U) + Pr*Prt*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(-2*(dUdx - 2*dVdy)*(dNTLdy*RHO + dRHOdy*NTL)*V + 2*(2*dUdx - dVdy)*(dNTLdx*RHO + dRHOdx*NTL)*U + 3*(dUdy + dVdx)*(dNTLdx*RHO + dRHOdx*NTL)*V + 3*(dUdy + dVdx)*(dNTLdy*RHO + dRHOdy*NTL)*U)*pow(NTL, 2)*pow(RHO, 2) + Gamma*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*(-d2Edxx + 1.0*d2Udxx*U + 1.0*d2Vdxx*V + 1.0*pow(dUdx, 2) + 1.0*pow(dVdx, 2)) + (Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*(-d2Edyy + 1.0*d2Udyy*U + 1.0*d2Vdyy*V + 1.0*pow(dUdy, 2) + 1.0*pow(dVdy, 2)) + (-dEdx + 1.0*dUdx*U + 1.0*dVdx*V)*(4*Pr*dNTLdx*NTL*pow(RHO, 2) + 4*Pr*dRHOdx*pow(NTL, 2)*RHO + 3*Prt*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2) + (-dEdy + 1.0*dUdy*U + 1.0*dVdy*V)*(4*Pr*dNTLdy*NTL*pow(RHO, 2) + 4*Pr*dRHOdy*pow(NTL, 2)*RHO + 3*Prt*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) - 3*Gamma*((dNTLdx*RHO + dRHOdx*NTL)*(-dEdx + 1.0*dUdx*U + 1.0*dVdx*V) + (dNTLdy*RHO + dRHOdy*NTL)*(-dEdy + 1.0*dUdy*U + 1.0*dVdy*V))*(Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*pow(NTL, 2)*pow(RHO, 2))/(Pr*Prt*pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2));
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;
      
    // qnt_id=4: NTL
    case 4:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          NTL = 2.71828182845905*NTL_max*pow(sig_nu, 2)*pow(y, 2)*exp(-pow(sig_nu, 2)*pow(y, 2)/pow(x, 2))/pow(x, 2);
          qnt = NTL;
          break;
          
        // qnt_type=1: dQdx
        case 1:                        
          dNTLdx = 5.43656365691809*NTL_max*pow(sig_nu, 4)*pow(y, 4)*exp(-pow(sig_nu, 2)*pow(y, 2)/pow(x, 2))/pow(x, 5) - 5.43656365691809*NTL_max*pow(sig_nu, 2)*pow(y, 2)*exp(-pow(sig_nu, 2)*pow(y, 2)/pow(x, 2))/pow(x, 3);
          qnt  = dNTLdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dNTLdy = -5.43656365691809*NTL_max*pow(sig_nu, 4)*pow(y, 3)*exp(-pow(sig_nu, 2)*pow(y, 2)/pow(x, 2))/pow(x, 4) + 5.43656365691809*NTL_max*pow(sig_nu, 2)*y*exp(-pow(sig_nu, 2)*pow(y, 2)/pow(x, 2))/pow(x, 2);
          qnt  = dNTLdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2NTLdxx = NTL_max*pow(sig_nu, 2)*pow(y, 2)*(10.8731273138362*pow(sig_nu, 4)*pow(y, 4)/pow(x, 4) - 38.0559455984266*pow(sig_nu, 2)*pow(y, 2)/pow(x, 2) + 16.3096909707543)*exp(-pow(sig_nu, 2)*pow(y, 2)/pow(x, 2))/pow(x, 4); 
          qnt    = d2NTLdxx;
          break;      

        // qnt_type=4: d2Qdxy
//         case 4:    
//           d2NTLdxy = 4*sig_v*y*(pow(sig_v, 2)*pow(y, 2)/pow(x, 2) - 1)*exp(-pow(sig_v, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 3)); 
//           qnt    = d2NTLdxy;
//           break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2NTLdyy = NTL_max*pow(sig_nu, 2)*(10.8731273138362*pow(sig_nu, 4)*pow(y, 4)/pow(x, 4) - 27.1828182845905*pow(sig_nu, 2)*pow(y, 2)/pow(x, 2) + 5.43656365691809)*exp(-pow(sig_nu, 2)*pow(y, 2)/pow(x, 2))/pow(x, 2);  
          qnt    = d2NTLdyy;
          break;                

     // qnt_type=10: forcing function
        case 10:

          FFQ_NTL_INV          = MS_qnt_2D(MS_no,qt_FFQ_inv,    qi_SA,x,y);
          FFQ_NTL_VIS          = MS_qnt_2D(MS_no,qt_FFQ_vis,    qi_SA,x,y);
          FFQ_NTL_SRC_PRODUCT  = MS_qnt_2D(MS_no,qt_FFQ_SA_prod,qi_SA,x,y);
          FFQ_NTL_SRC_DESTRUCT = MS_qnt_2D(MS_no,qt_FFQ_SA_dest,qi_SA,x,y);
          FFQ_NTL_SRC_DISTRIB  = MS_qnt_2D(MS_no,qt_FFQ_SA_dist,qi_SA,x,y);
          FFQ_NTL_SRC_CONSERV  = MS_qnt_2D(MS_no,qt_FFQ_SA_cons,qi_SA,x,y);
          
          FFQ_NTL = FFQ_NTL_INV + FFQ_NTL_VIS + FFQ_NTL_SRC_PRODUCT  + FFQ_NTL_SRC_DESTRUCT + FFQ_NTL_SRC_DISTRIB + FFQ_NTL_SRC_CONSERV; 
          qnt    = FFQ_NTL;
          break;
       
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U        = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          
          V        = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                       

          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                

          FFQ_NTL_INV = dNTLdx*RHO*U + dNTLdy*RHO*V + dRHOdx*NTL*U + dRHOdy*NTL*V + dUdx*NTL*RHO + dVdy*NTL*RHO; 
          qnt    = FFQ_NTL_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
        
          NTL      = MS_qnt_2D(MS_no,qt_Q,     qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Ntl,x,y);                       
          d2NTLdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Ntl,x,y);
          d2NTLdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Ntl,x,y);
                    
          FFQ_NTL_VIS = -(d2NTLdxx*(mu + NTL*RHO) + d2NTLdyy*(mu + NTL*RHO) + dNTLdx*(dNTLdx*RHO + dRHOdx*NTL) + dNTLdy*(dNTLdy*RHO + dRHOdy*NTL))/sigma;           
          qnt = FFQ_NTL_VIS;
          break;
          
         // qnt_type=13: Production term in the SA source term
        case 13:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);             

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);          
          
          if (d>0.){                                      
            // The negative portion of the SA model is activated based on two criteria:
            // for the S_tilda (modified vorticity used in the P and D terms), is computed differently
            // based on "if(Sbar>=-cv2*S)".
            // As for all the other values/terms/quantities, the negative portion
            // is activated directly upon the sign of the NTL (SA Working variable)
            S_MS    = fabs(dVdx-dUdy);
            nu_MS   = mu/RHO;
            chi_MS  = NTL/nu_MS;
            fv1_MS  = pow(chi_MS,3)/(pow(chi_MS,3)+pow(c_v1,3));
            fv2_MS  = 1.-(chi_MS/(1.+chi_MS*fv1_MS));
            Sbar_MS = fv2_MS*NTL/((kappa*kappa)*(d*d)); 
            
            if(Sbar_MS>=-c_v2*S_MS){
              if(NTL>=0.){            
                FFQ_NTL_SRC_PRODUCT = -c_b1*(-c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2)) + 1)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*NTL*RHO; 
              }
              else{
                printf("MMS for wall-bounded negative SA is not implemented yet\n");
                exit(1);
              }                
            }
            else
            {
              if(NTL>=0.){
                FFQ_NTL_SRC_PRODUCT = -c_b1*(-c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2)) + 1)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx))*NTL*RHO; 
              }
              else{
                printf("MMS for wall-bounded negative SA is not implemented yet\n");
                exit(1);
              }              
            }                      
                        
          }                            
          else
          {
           tau_w = dUdy*mu; 
            FFQ_NTL_SRC_PRODUCT = c_b1*c_t3*tau_w - c_b1*tau_w; 
          }
          
          qnt = FFQ_NTL_SRC_PRODUCT;
          break; 

        // qnt_type=14: Destruction term in the SA source term
        case 14:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);              

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);           
          
          if (d>0.){                                      
            // The negative portion of the SA model is activated based on two criteria:
            // for the S_tilda (modified vorticity used in the P and D terms), is computed differently
            // based on "if(Sbar>=-cv2*S)".
            // As for all the other values/terms/quantities, the negative portion
            // is activated directly upon the sign of the NTL (SA Working variable)
            S_MS    = fabs(dVdx-dUdy);
            nu_MS   = mu/RHO;
            chi_MS  = NTL/nu_MS;
            fv1_MS  = pow(chi_MS,3)/(pow(chi_MS,3)+pow(c_v1,3));
            fv2_MS  = 1.-(chi_MS/(1.+chi_MS*fv1_MS));
            Sbar_MS = fv2_MS*NTL/((kappa*kappa)*(d*d)); 
            
            if(Sbar_MS>=-c_v2*S_MS){
              if(NTL>=0.){            
                FFQ_NTL_SRC_DESTRUCT = (-c_b1*c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2))/pow(kappa, 2) + c_w1*pow((pow(c_w3, 6) + 1)/(pow(c_w3, 6) + pow(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6)), 1.0L/6.0L)*(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))))*pow(NTL, 2)*RHO/pow(d, 2);
              }
              else{
                printf("MMS for wall-bounded negative SA is not implemented yet\n");
                exit(1);
              }                
            }
            else
            {
              if(NTL>=0.){
                FFQ_NTL_SRC_DESTRUCT = (-c_b1*c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2))/pow(kappa, 2) + c_w1*pow((pow(c_w3, 6) + 1)/(pow(c_w3, 6) + pow(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx)))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx)))), 6)), 1.0L/6.0L)*(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx)))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx))))))*pow(NTL, 2)*RHO/pow(d, 2); 
              }
              else{
                printf("MMS for wall-bounded negative SA is not implemented yet\n");
                exit(1);
              }              
            }                      
                        
          }                            
          else
          {
            tau_w = dUdy*mu; 
            FFQ_NTL_SRC_DESTRUCT = -c_b1*c_t3*tau_w + c_w1*pow(kappa, 2)*tau_w;
          } 
          qnt = FFQ_NTL_SRC_DESTRUCT;
          break;  

        // qnt_type=15: Distribtuion term in the SA source term 
        case 15:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                 
          
          FFQ_NTL_SRC_DISTRIB = -c_b2*(pow(dNTLdx, 2) + pow(dNTLdy, 2))*RHO/sigma;           
          qnt = FFQ_NTL_SRC_DISTRIB;
          break;            
          
        // qnt_type=13: Conservation term in the SA source term 
        case 16:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
        
          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                       
                    
          FFQ_NTL_SRC_CONSERV = (dNTLdx*dRHOdx + dNTLdy*dRHOdy)*(mu/RHO + NTL)/sigma;
          qnt = FFQ_NTL_SRC_CONSERV;
          break;               
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;       
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
       case 0:              
          P = p_0 + 0.5*log(-pow(x, 2) + 2*x + 0.25)*log(4*pow(y, 3) - 3*pow(y, 2) + 1.25); 
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
	  dPdx = 0.5*(-2*x + 2)*log(4*pow(y, 3) - 3*pow(y, 2) + 1.25)/(-pow(x, 2) + 2*x + 0.25); 
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = 0.5*(12*pow(y, 2) - 6*y)*log(-pow(x, 2) + 2*x + 0.25)/(4*pow(y, 3) - 3*pow(y, 2) + 1.25); 
          qnt  = dPdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:    
          d2Pdxx = -(2.0*pow(x - 1, 2)/(-pow(x, 2) + 2*x + 0.25) + 1.0)*log(4*pow(y, 3) - 3*pow(y, 2) + 1.25)/(-pow(x, 2) + 2*x + 0.25); 
          qnt = d2Pdxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:                        
          d2Pdyy = (-18.0*pow(y, 2)*pow(2*y - 1, 2)/(4*pow(y, 3) - 3*pow(y, 2) + 1.25) + 12.0*y - 3.0)*log(-pow(x, 2) + 2*x + 0.25)/(4*pow(y, 3) - 3*pow(y, 2) + 1.25);
          qnt    = d2Pdyy;
          break;               

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;  
      
    // qnt_id=10: d
    case 10:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          d   = y;
          qnt = d;
          break;
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;        
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}

double MS_9(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution and forcing functions      
      
    // MS-9: Realistic RANS-SA Verification case, Oliver et al. (2012)
      
*/  

  
// Do not change these values of model parameters, instead set the values in the solver to these.
// This is important for correct verification as well as cross-solvers comparability.  
  const double Pr     = 0.7    ;
  const double mu     = 1e-4   ; // 2.67861904719577e-6; // for the non-dimensional version
  const double Gamma  = 1.4    ;  

  const double C_cf   = 0.0270 ;
  const double eta_1  = 11.0   ;
  const double b      = 0.330  ;
  const double C      = 5.0    ;
  const double eta_v  = 30.0   ; 
  const double T_inf  = 250.0  ; // 1.0    ;  // for the non-dimensional version
  const double M_inf  = 0.80   ;
  const double r_T    = 0.90   ;
  const double alpha  = 5.0    ; // 1.86663348236639e-2;  // for the non-dimensional version
   
  const double p0     = 1.0e4  ; // 1.0    ;  // for the non-dimensional version
  const double R      = 287.0  ; // 1.0    ;  // for the non-dimensional version
    
  double u_inf  = M_inf * sqrt(Gamma * R * T_inf);
  double T_w    = T_inf *(1. + r_T * (Gamma-1.)/2. * pow(M_inf,2.));
  double ro_w   = p0/(R * T_w);
  double ro_inf = p0/(R * T_inf);
  double v_w    = mu/ro_w;
  double A      = sqrt(1. - T_inf/T_w);
  double F_c    = (T_w/T_inf - 1.)/pow(asin(A),2.);

  // RANS-SA variables
  const double sigma    = 2./3.;
  const double Prt      = 0.9  ;  
  const double c_b1     = 0.1355;
  const double c_b2     = 0.622;
  const double c_t3     = 1.2;
  const double c_t4     = 0.5;       
  const double c_v1     = 7.1;
  const double c_v2     = 0.7;
  const double c_v3     = 0.9;
  const double kappa    = 0.41;     
  const double c_w1     = c_b1/(kappa*kappa)+(1.0+c_b2)/sigma;
  const double c_w2     = 0.3;
  const double c_w3     = 2.0;  
  const double r_lim    = 10.0;
  
  double S_MS,Sbar_MS,nu_MS,chi_MS,fv1_MS,fv2_MS;    
  double RHO,U,V,E,P,NTL;
  double d;
  double dRHOdx,dRHOdy,d2RHOdxx,d2RHOdyy;  
  double dUdx,dUdy,d2Udxx,d2Udxy,d2Udyy;
  double dVdx,dVdy,d2Vdxx,d2Vdxy,d2Vdyy;  
  double dEdx,dEdy,d2Edxx,d2Edyy; // double d2Edxy;
  double dNTLdx,dNTLdy,d2NTLdxx,d2NTLdyy; // double d2NTLdxy; 
  double dPdx,dPdy,d2Pdxx,d2Pdyy;
  double FFQ_RHO_INV,FFQ_U_INV,FFQ_V_INV,FFQ_E_INV,FFQ_NTL_INV;
  double FFQ_RHO_VIS,FFQ_U_VIS,FFQ_V_VIS,FFQ_E_VIS,FFQ_NTL_VIS;  
  double FFQ_RHO,FFQ_U,FFQ_V,FFQ_E,FFQ_NTL;
  double FFQ_NTL_SRC_PRODUCT,FFQ_NTL_SRC_DESTRUCT,FFQ_NTL_SRC_DISTRIB,FFQ_NTL_SRC_CONSERV;   
  double tau_w;  
  
  double qnt;
  int MS_no;
  
  MS_no = 9;
  
  switch(qnt_id)
  {
       
    // qnt_id=0: ro
    case 0:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO = p0/(R*T_inf*((1.0L/2.0L)*pow(M_inf, 2)*r_T*(1 - pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/pow(A, 2))*(Gamma - 1) + 1));  
          qnt = RHO;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dRHOdx = pow(M_inf, 2)*p0*r_T*(Gamma - 1)*((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/28.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(v_w*x*((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)) + (C - log(kappa)/kappa)*((1.0L/28.0L)*b*pow(u_inf, 2)*pow(y, 2)*(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*pow(v_w, 2)*x) + (1.0L/28.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/(eta_1*v_w*x) - 1.0L/28.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w*x))) - 1.0L/28.0L*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)/x)*sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))/(pow(A, 2)*R*T_inf*pow((1.0L/2.0L)*pow(M_inf, 2)*r_T*(1 - pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/pow(A, 2))*(Gamma - 1) + 1, 2)); 
          qnt    = dRHOdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dRHOdy = (1.0L/2.0L)*sqrt(2)*pow(M_inf, 2)*p0*r_T*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*(Gamma - 1)*(-1.0L/2.0L*sqrt(2)*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(v_w*((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)) + (C - log(kappa)/kappa)*(-1.0L/2.0L*b*pow(u_inf, 2)*y*C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*pow(v_w, 2)) - 1.0L/2.0L*sqrt(2)*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/(eta_1*v_w) + (1.0L/2.0L)*sqrt(2)*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)))*sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))/(A*R*T_inf*pow((1.0L/2.0L)*pow(M_inf, 2)*r_T*(1 - pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/pow(A, 2))*(Gamma - 1) + 1, 2));
          qnt    = dRHOdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:
          d2RHOdxx = (1.0L/392.0L)*pow(M_inf, 2)*p0*r_T*(Gamma - 1)*(-C_cf*pow(u_inf*y*(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w - (C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) + 2*log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa, 2)*pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + C_cf*pow(u_inf*y*(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w - (C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) + 2*log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa, 2)*pow(cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*(-4*u_inf*y*(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w + u_inf*y*(8*C_cf*kappa*u_inf*y/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)*pow(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0, 2)) - 60*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(sqrt(2)*C_cf*pow(b, 2)*pow(u_inf, 2)*pow(y, 2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*pow(v_w, 2)*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) - 34*C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + 2*C_cf*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/(F_c*eta_1*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) - 30*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + 30*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w + 30*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - 60*log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)*sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))/A + 4*C_cf*pow(M_inf, 2)*r_T*(Gamma - 1)*pow(u_inf*y*(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w - (C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) + 2*log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa, 2)*pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)*pow(cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/(pow(A, 2)*F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)*(pow(M_inf, 2)*r_T*(1 - pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/pow(A, 2))*(Gamma - 1) + 2)))/(R*T_inf*pow(x, 2)*pow(pow(M_inf, 2)*r_T*(1 - pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/pow(A, 2))*(Gamma - 1) + 2, 2)); 
          qnt    = d2RHOdxx;
          break;      

        // qnt_type=4: d2Qdxy
//        case 4:    
//          qnt    = d2RHOdxy;
//          break;
          
        // qnt_type=5: d2Qdyy
        case 5:
          d2RHOdyy = C_cf*pow(M_inf, 2)*p0*r_T*pow(u_inf, 2)*(Gamma - 1)*(-1.0L/2.0L*pow(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1, 2)*pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2) + (1.0L/2.0L)*pow(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1, 2)*pow(cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2) + (1.0L/2.0L)*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*(8*kappa/pow(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0, 2) + (C - log(kappa)/kappa)*(sqrt(2)*pow(b, 2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/v_w - 4*b*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w) + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/eta_1)/eta_1)*sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))/A + 2*pow(M_inf, 2)*r_T*(Gamma - 1)*pow(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1, 2)*pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)*pow(cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/(pow(A, 2)*(pow(M_inf, 2)*r_T*(1 - pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/pow(A, 2))*(Gamma - 1) + 2)))/(F_c*R*T_inf*pow(v_w, 2)*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)*pow(pow(M_inf, 2)*r_T*(1 - pow(sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)), 2)/pow(A, 2))*(Gamma - 1) + 2, 2));
          qnt    = d2RHOdyy;
          break;           

        // qnt_type=10: forcing function
        case 10:

          FFQ_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          FFQ_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          FFQ_RHO = FFQ_RHO_INV + FFQ_RHO_VIS; 
          qnt    = FFQ_RHO;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);          
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);
          
          FFQ_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = FFQ_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          FFQ_RHO_VIS = 0.;
          qnt    = FFQ_RHO_VIS;
          break; 
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  

      }
      break;

    // qnt_id=1: U
    case 1:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          U = -u_inf*sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))/A; 
          qnt = U;
          break;
          
        // qnt_type=1: dQdx
        case 1:    
          dUdx = -u_inf*((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/28.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(v_w*x*((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)) + (C - log(kappa)/kappa)*((1.0L/28.0L)*b*pow(u_inf, 2)*pow(y, 2)*(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*pow(v_w, 2)*x) + (1.0L/28.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/(eta_1*v_w*x) - 1.0L/28.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w*x))) - 1.0L/28.0L*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)/x)*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))/A;
          qnt  = dUdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dUdy = -1.0L/2.0L*sqrt(2)*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*(-1.0L/2.0L*sqrt(2)*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(v_w*((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)) + (C - log(kappa)/kappa)*(-1.0L/2.0L*b*pow(u_inf, 2)*y*C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*pow(v_w, 2)) - 1.0L/2.0L*sqrt(2)*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/(eta_1*v_w) + (1.0L/2.0L)*sqrt(2)*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)))*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((C - log(kappa)/kappa)*(-1 + exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + (1.0L/2.0L)*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa));
          qnt  = dUdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2Udxx = (1.0L/1568.0L)*u_inf*(A*C_cf*pow(u_inf*y*(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w - (C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) + 2*log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa, 2)*sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*(-4*u_inf*y*(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w + u_inf*y*(8*C_cf*kappa*u_inf*y/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)*pow(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0, 2)) - 60*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(sqrt(2)*C_cf*pow(b, 2)*pow(u_inf, 2)*pow(y, 2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*pow(v_w, 2)*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) - 34*C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + 2*C_cf*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/(F_c*eta_1*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) - 30*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + 30*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w + 30*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - 60*log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)))/pow(x, 2); 
          qnt    = d2Udxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:    
          d2Udxy = -1.0L/112.0L*pow(u_inf, 2)*(A*C_cf*(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)*(u_inf*y*(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)/v_w - (C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) + 2*log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)*sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa))/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*(-8*C_cf*kappa*u_inf*y/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)*pow(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0, 2)) + 8*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + 2*(C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1 + (C - log(kappa)/kappa)*(-sqrt(2)*C_cf*pow(b, 2)*pow(u_inf, 2)*pow(y, 2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*pow(v_w, 2)*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + 6*C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) - 2*C_cf*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/(F_c*eta_1*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + 2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - 2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1)*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)))/(v_w*x); 
          qnt    = d2Udxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:    
          d2Udyy = (1.0L/8.0L)*C_cf*pow(u_inf, 3)*(A*pow(2*sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0) + (C - log(kappa)/kappa)*(C_cf*b*u_inf*y*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(F_c*v_w*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)) + sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w))/eta_1, 2)*sin((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)) - sqrt(2)*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*(8*kappa/pow(sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 2.0, 2) + (C - log(kappa)/kappa)*(sqrt(2)*pow(b, 2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/v_w - 4*b*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w) + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w))/eta_1)/eta_1)*cos((1.0L/2.0L)*sqrt(2)*A*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*((1.0L/2.0L)*(C - log(kappa)/kappa)*(-2 + 2*exp(-1.0L/2.0L*sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/(eta_1*v_w)) + sqrt(2)*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))*exp(-1.0L/2.0L*sqrt(2)*b*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w)/(eta_1*v_w)) - log((1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/v_w + 1.0)/kappa)))/(F_c*pow(v_w, 2)*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)); 
          qnt    = d2Udyy;
          break;                

        // qnt_type=10: forcing function
        case 10:

          FFQ_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          FFQ_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          FFQ_U = FFQ_U_INV + FFQ_U_VIS; 
          qnt    = FFQ_U;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);      

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                      
                    
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);                    

          FFQ_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U; 
          qnt    = FFQ_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);          
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);          

          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);  
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);            
                 
          FFQ_U_VIS = (-2.0L/3.0L*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((2*d2Udxx - d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (2*dUdx - dVdy)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)) - (pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udyy + d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdy + dVdx)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(2*(2*dUdx - dVdy)*(dNTLdx*RHO + dRHOdx*NTL) + 3*(dUdy + dVdx)*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2))/pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2);          
          qnt    = FFQ_U_VIS;
          break;                          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;           
      }
      break;
      
    // qnt_id=2: V
    case 2:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          V =  (1.0L/28.0L)*sqrt(2)*eta_v*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/x; 
          qnt = V;
          break;
          
        // qnt_type=1: dQdx
        case 1:              
          dVdx = -15.0L/392.0L*sqrt(2)*eta_v*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/pow(x, 2);
          qnt  = dVdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dVdy = (1.0L/28.0L)*sqrt(2)*eta_v*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/x;
          qnt  = dVdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2Vdxx = (435.0L/5488.0L)*sqrt(2)*eta_v*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/pow(x, 3);
          qnt    = d2Vdxx;
          break;      

        // qnt_type=4: d2Qdxy
        case 4:              
          d2Vdxy = -15.0L/392.0L*sqrt(2)*eta_v*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/pow(x, 2);
          qnt    = d2Vdxy;
          break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2Vdyy = 0; 
          qnt    = d2Vdyy;
          break;                

      // qnt_type=10: forcing function
        case 10:

          FFQ_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          FFQ_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          FFQ_V = FFQ_V_INV + FFQ_V_VIS; 
          qnt    = FFQ_V;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);             
                  
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);             
         
          FFQ_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V; 
          qnt    = FFQ_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);          

          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);   
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);                          
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y); 
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);
          
          FFQ_V_VIS = (-(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udxy + d2Vdxx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdy + dVdx)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (2.0L/3.0L)*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udxy - 2*d2Vdyy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdx - 2*dVdy)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(-2*(dUdx - 2*dVdy)*(dNTLdy*RHO + dRHOdy*NTL) + 3*(dUdy + dVdx)*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2))/pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2);
          
          qnt    = FFQ_V_VIS;
          break;          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;  
      
    // qnt_id=3: E
    case 3:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,x,y);
          P      = MS_qnt_2D(MS_no,qt_Q,qi_Prsr,x,y);
          
          E = 0.5*pow(U, 2.0) + 0.5*pow(V, 2.0) + P/((Gamma - 1.0)*RHO); 
          qnt = E;
          break;
          
        // qnt_type=1: dQdx
        case 1:  
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);              
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);          
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);
                    
          dEdx = 1.0*U*dUdx + 1.0*V*dVdx - P*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdx/((Gamma - 1.0)*RHO); 
          qnt = dEdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);    

          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);                       
          
          dEdy = 1.0*U*dUdy + 1.0*V*dVdy - P*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)) + dPdy/((Gamma - 1.0)*RHO); 
          qnt = dEdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:  
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Rho,x,y);
          d2RHOdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Rho,x,y);
          
          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          d2Udxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y); 
          d2Vdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);            
                   
          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);
          dPdx     = MS_qnt_2D(MS_no,qt_dQdx,  qi_Prsr,x,y);
          d2Pdxx   = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Prsr,x,y); 
                    
          d2Edxx = 1.0*U*d2Udxx + 1.0*V*d2Vdxx + 1.0*pow(dUdx, 2.0) + 1.0*pow(dVdx, 2.0) - P*d2RHOdxx/((Gamma - 1.0)*pow(RHO, 2.0)) + 2.0*P*pow(dRHOdx, 2.0)/((Gamma - 1.0)*pow(RHO, 3.0)) + d2Pdxx/((Gamma - 1.0)*RHO) - 2.0*dPdx*dRHOdx/((Gamma - 1.0)*pow(RHO, 2.0));
          qnt = d2Edxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:
          RHO      = MS_qnt_2D(MS_no,qt_Q,     qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Rho,x,y);
          d2RHOdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Rho,x,y);          

          U        = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);
                   
          V        = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          P        = MS_qnt_2D(MS_no,qt_Q,     qi_Prsr,x,y);      
          dPdy     = MS_qnt_2D(MS_no,qt_dQdy,  qi_Prsr,x,y);  
          d2Pdyy   = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Prsr,x,y);                       
                   
          d2Edyy = 1.0*U*d2Udyy + 1.0*V*d2Vdyy + 1.0*pow(dUdy, 2.0) + 1.0*pow(dVdy, 2.0) - P*d2RHOdyy/((Gamma - 1.0)*pow(RHO, 2.0)) + 2.0*P*pow(dRHOdy, 2.0)/((Gamma - 1.0)*pow(RHO, 3.0)) + d2Pdyy/((Gamma - 1.0)*RHO) - 2.0*dPdy*dRHOdy/((Gamma - 1.0)*pow(RHO, 2.0)); 
          qnt    = d2Edyy;
          break;        
          
     // qnt_type=10: forcing function
        case 10:

          FFQ_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          FFQ_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          FFQ_E = FFQ_E_INV + FFQ_E_VIS; 
          qnt    = FFQ_E;
          break;
          
        // qnt_type=11: inviscid component of the forcing function 
        case 11:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                            
          
          P      = MS_qnt_2D(MS_no,qt_Q,   qi_Prsr,x,y);
          dPdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Prsr,x,y);          
          dPdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Prsr,x,y);           
          
          E      = MS_qnt_2D(MS_no,qt_Q,   qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_E,x,y);       

          FFQ_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;
          qnt    = FFQ_E_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:

          RHO    = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U      = MS_qnt_2D(MS_no,qt_Q,     qi_U,x,y);
          dUdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_U,x,y);
          dUdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_U,x,y);
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);

          V      = MS_qnt_2D(MS_no,qt_Q,     qi_V,x,y);
          dVdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_V,x,y);
          dVdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_V,x,y);    
          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);             

          E      = MS_qnt_2D(MS_no,qt_Q,     qi_E,x,y);
          dEdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_E,x,y);          
          dEdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_E,x,y);
          d2Edxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_E,x,y);
          d2Edyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_E,x,y);
          
          NTL    = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y); 
          dNTLdy = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                      
          
          
          FFQ_E_VIS = ((2.0L/3.0L)*Pr*Prt*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(-dUdx*(2*dUdx - dVdy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + dVdy*(dUdx - 2*dVdy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) - (2*d2Udxx - d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*U + (d2Udxy - 2*d2Vdyy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*V + (dUdx - 2*dVdy)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)*V - (2*dUdx - dVdy)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)*U) - Pr*Prt*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(dUdy*(dUdy + dVdx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + dVdx*(dUdy + dVdx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (d2Udxy + d2Vdxx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*V + (d2Udyy + d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*U + (dUdy + dVdx)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)*V + (dUdy + dVdx)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)*U) + Pr*Prt*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(-2*(dUdx - 2*dVdy)*(dNTLdy*RHO + dRHOdy*NTL)*V + 2*(2*dUdx - dVdy)*(dNTLdx*RHO + dRHOdx*NTL)*U + 3*(dUdy + dVdx)*(dNTLdx*RHO + dRHOdx*NTL)*V + 3*(dUdy + dVdx)*(dNTLdy*RHO + dRHOdy*NTL)*U)*pow(NTL, 2)*pow(RHO, 2) + Gamma*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*(-d2Edxx + 1.0*d2Udxx*U + 1.0*d2Vdxx*V + 1.0*pow(dUdx, 2) + 1.0*pow(dVdx, 2)) + (Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*(-d2Edyy + 1.0*d2Udyy*U + 1.0*d2Vdyy*V + 1.0*pow(dUdy, 2) + 1.0*pow(dVdy, 2)) + (-dEdx + 1.0*dUdx*U + 1.0*dVdx*V)*(4*Pr*dNTLdx*NTL*pow(RHO, 2) + 4*Pr*dRHOdx*pow(NTL, 2)*RHO + 3*Prt*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2) + (-dEdy + 1.0*dUdy*U + 1.0*dVdy*V)*(4*Pr*dNTLdy*NTL*pow(RHO, 2) + 4*Pr*dRHOdy*pow(NTL, 2)*RHO + 3*Prt*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) - 3*Gamma*((dNTLdx*RHO + dRHOdx*NTL)*(-dEdx + 1.0*dUdx*U + 1.0*dVdx*V) + (dNTLdy*RHO + dRHOdy*NTL)*(-dEdy + 1.0*dUdy*U + 1.0*dVdy*V))*(Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*pow(NTL, 2)*pow(RHO, 2))/(Pr*Prt*pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2)); 
          qnt    = FFQ_E_VIS;
          break;             

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;  
          
      }
      break;
      
    // qnt_id=4: NTL
    case 4:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:          
          NTL = -alpha*pow(y, 2) + (1.0L/2.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)));
          qnt = NTL;
          break;
          
        // qnt_type=1: dQdx
        case 1:                        
          dNTLdx = -1.0L/28.0L*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/x;
          qnt  = dNTLdx;
          break;          

        // qnt_type=2: dQdy
        case 2:    
          dNTLdy = -2*alpha*y + (1.0L/2.0L)*sqrt(2)*kappa*u_inf*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)));
          qnt  = dNTLdy;
          break;      

        // qnt_type=3: d2Qdxx
        case 3:    
          d2NTLdxx = (15.0L/392.0L)*sqrt(2)*kappa*u_inf*y*sqrt(C_cf/(F_c*pow(ro_inf*u_inf*x/(F_c*mu), 1.0L/7.0L)))/pow(x, 2);
          qnt    = d2NTLdxx;
          break;      

        // qnt_type=4: d2Qdxy
//         case 4:    
//           d2NTLdxy = 4*sig_v*y*(pow(sig_v, 2)*pow(y, 2)/pow(x, 2) - 1)*exp(-pow(sig_v, 2)*pow(y, 2)/pow(x, 2))/(sqrt(piMS)*pow(x, 3)); 
//           qnt    = d2NTLdxy;
//           break;
          
        // qnt_type=5: d2Qdyy
        case 5:              
          d2NTLdyy = -2*alpha;
          qnt    = d2NTLdyy;
          break;                

     // qnt_type=10: forcing function
        case 10:

          FFQ_NTL_INV          = MS_qnt_2D(MS_no,qt_FFQ_inv,    qi_SA,x,y);
          FFQ_NTL_VIS          = MS_qnt_2D(MS_no,qt_FFQ_vis,    qi_SA,x,y);
          FFQ_NTL_SRC_PRODUCT  = MS_qnt_2D(MS_no,qt_FFQ_SA_prod,qi_SA,x,y);
          FFQ_NTL_SRC_DESTRUCT = MS_qnt_2D(MS_no,qt_FFQ_SA_dest,qi_SA,x,y);
          FFQ_NTL_SRC_DISTRIB  = MS_qnt_2D(MS_no,qt_FFQ_SA_dist,qi_SA,x,y);
          FFQ_NTL_SRC_CONSERV  = MS_qnt_2D(MS_no,qt_FFQ_SA_cons,qi_SA,x,y);
          
          FFQ_NTL = FFQ_NTL_INV + FFQ_NTL_VIS + FFQ_NTL_SRC_PRODUCT  + FFQ_NTL_SRC_DESTRUCT + FFQ_NTL_SRC_DISTRIB + FFQ_NTL_SRC_CONSERV; 
          qnt    = FFQ_NTL;
          break;
       
        // qnt_type=11: inviscid component of the forcing function 
        case 11:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);

          U        = MS_qnt_2D(MS_no,qt_Q,   qi_U,x,y);
          dUdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_U,x,y);
          
          V        = MS_qnt_2D(MS_no,qt_Q,   qi_V,x,y);
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);                       

          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                

          FFQ_NTL_INV = dNTLdx*RHO*U + dNTLdy*RHO*V + dRHOdx*NTL*U + dRHOdy*NTL*V + dUdx*NTL*RHO + dVdy*NTL*RHO; 
          qnt    = FFQ_NTL_INV;
          break;  
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
        
          NTL      = MS_qnt_2D(MS_no,qt_Q,     qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,  qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,  qi_Ntl,x,y);                       
          d2NTLdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_Ntl,x,y);
          d2NTLdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_Ntl,x,y);
                    
          FFQ_NTL_VIS = -(d2NTLdxx*(mu + NTL*RHO) + d2NTLdyy*(mu + NTL*RHO) + dNTLdx*(dNTLdx*RHO + dRHOdx*NTL) + dNTLdy*(dNTLdy*RHO + dRHOdy*NTL))/sigma;         
          qnt = FFQ_NTL_VIS;
          break;
        // qnt_type=13: Production term in the SA source term
        case 13:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);             

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);          
          
          if (d>0.){                                      
            // The negative portion of the SA model is activated based on two criteria:
            // for the S_tilda (modified vorticity used in the P and D terms), is computed differently
            // based on "if(Sbar>=-cv2*S)".
            // As for all the other values/terms/quantities, the negative portion
            // is activated directly upon the sign of the NTL (SA Working variable)
            S_MS    = fabs(dVdx-dUdy);
            nu_MS   = mu/RHO;
            chi_MS  = NTL/nu_MS;
            fv1_MS  = pow(chi_MS,3)/(pow(chi_MS,3)+pow(c_v1,3));
            fv2_MS  = 1.-(chi_MS/(1.+chi_MS*fv1_MS));
            Sbar_MS = fv2_MS*NTL/((kappa*kappa)*(d*d)); 
            
            if(Sbar_MS>=-c_v2*S_MS){
              if(NTL>=0.){            
                FFQ_NTL_SRC_PRODUCT = -c_b1*(-c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2)) + 1)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*NTL*RHO;  
              }
              else{
                printf("MMS for wall-bounded negative SA is not implemented yet\n");
                exit(1);
              }                
            }
            else
            {
              if(NTL>=0.){
                FFQ_NTL_SRC_PRODUCT = -c_b1*(-c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2)) + 1)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx))*NTL*RHO;
              }
              else{
                printf("MMS for wall-bounded negative SA is not implemented yet\n");
                exit(1);
              }              
            }                      
                        
          }                            
          else
          {
            tau_w = dUdy*mu; 
            FFQ_NTL_SRC_PRODUCT = c_b1*c_t3*tau_w - c_b1*tau_w; 
          }
          
          qnt = FFQ_NTL_SRC_PRODUCT;
          break; 

        // qnt_type=14: Destruction term in the SA source term
        case 14:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);
          
          dVdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_V,x,y);              

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);           
          
          if (d>0.){                                      
            // The negative portion of the SA model is activated based on two criteria:
            // for the S_tilda (modified vorticity used in the P and D terms), is computed differently
            // based on "if(Sbar>=-cv2*S)".
            // As for all the other values/terms/quantities, the negative portion
            // is activated directly upon the sign of the NTL (SA Working variable)
            S_MS    = fabs(dVdx-dUdy);
            nu_MS   = mu/RHO;
            chi_MS  = NTL/nu_MS;
            fv1_MS  = pow(chi_MS,3)/(pow(chi_MS,3)+pow(c_v1,3));
            fv2_MS  = 1.-(chi_MS/(1.+chi_MS*fv1_MS));
            Sbar_MS = fv2_MS*NTL/((kappa*kappa)*(d*d)); 
            
            if(Sbar_MS>=-c_v2*S_MS){
              if(NTL>=0.){            
                FFQ_NTL_SRC_DESTRUCT = (-c_b1*c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2))/pow(kappa, 2) + c_w1*pow((pow(c_w3, 6) + 1)/(pow(c_w3, 6) + pow(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6)), 1.0L/6.0L)*(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))))*pow(NTL, 2)*RHO/pow(d, 2);
              }
              else{
                printf("MMS for wall-bounded negative SA is not implemented yet\n");
                exit(1);
              }                
            }
            else
            {
              if(NTL>=0.){
                FFQ_NTL_SRC_DESTRUCT = (-c_b1*c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2))/pow(kappa, 2) + c_w1*pow((pow(c_w3, 6) + 1)/(pow(c_w3, 6) + pow(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx)))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx)))), 6)), 1.0L/6.0L)*(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx)))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*((pow(c_v2, 2)*fabs(dUdy - dVdx) + c_v3*(1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*fabs(dUdy - dVdx)/((-2*c_v2 + c_v3)*fabs(dUdy - dVdx) - (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))) + fabs(dUdy - dVdx))))))*pow(NTL, 2)*RHO/pow(d, 2);
              }
              else{
                printf("MMS for wall-bounded negative SA is not implemented yet\n");
                exit(1);
              }              
            }                      
                        
          }                            
          else
          {
            tau_w = dUdy*mu; 
            FFQ_NTL_SRC_DESTRUCT = -c_b1*c_t3*tau_w + c_w1*pow(kappa, 2)*tau_w; 
          } 
          qnt = FFQ_NTL_SRC_DESTRUCT;
          break;  

        // qnt_type=15: Distribtuion term in the SA source term 
        case 15:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                 
          
          FFQ_NTL_SRC_DISTRIB = -c_b2*(pow(dNTLdx, 2) + pow(dNTLdy, 2))*RHO/sigma;
          qnt = FFQ_NTL_SRC_DISTRIB;
          break;            
          
        // qnt_type=13: Conservation term in the SA source term 
        case 16:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
        
          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                       
                    
          FFQ_NTL_SRC_CONSERV = (dNTLdx*dRHOdx + dNTLdy*dRHOdy)*(mu/RHO + NTL)/sigma; 
          qnt = FFQ_NTL_SRC_CONSERV;
          break;               
          
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);
          exit(1);
          break;          
      }
      break;       
 
    // qnt_id=9: P
    case 9:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
       case 0:              
          P = p0;
          qnt = P;
          break;
          
        // qnt_type=1: dQdx
        case 1:            
	  dPdx = 0.;
          qnt  = dPdx;
          break;          

        // qnt_type=2: dQdy
        case 2:          
          dPdy = 0.;
          qnt  = dPdy;
          break;  
          
        // qnt_type=3: d2Qdxx
        case 3:    
          d2Pdxx = 0.;
          qnt = d2Pdxx;
          break;      

//         // qnt_type=4: d2Qdxy
//         case 4:    
//           d2Edxy = 0.;
//           qnt = d2Edxy;
//           break;
//           
        // qnt_type=5: d2Qdyy
        case 5:                        
          d2Pdyy = 0.;
          qnt    = d2Pdyy;
          break;               

        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;  
      
    // qnt_id=10: d
    case 10:
      switch(qnt_type)
      {  
        // qnt_type=0: state variable Q
        case 0:    
          d   = y;
          qnt = d;
          break;
          
        default :
          printf("Invalid case in MS_qnt_2D\n");
          printf("qnt_type =%i\n",qnt_type);
          printf("qnt_id   =%i\n",qnt_id);          
          exit(1);
          break;            
      }
      break;        
      
    default :
      printf("Invalid case in MS_qnt_2D\n");
      printf("qnt_type =%i\n",qnt_type);
      printf("qnt_id   =%i\n",qnt_id);            
      exit(1);
      break;        
      
  }
  return(qnt);
}
