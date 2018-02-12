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

/* REFERENCE:
AIAA2017-3290                 
"Motivations and methods of verification for high-order RANS solvers and solutions"

Farshad Navah, McGill University; Sivakumaran Nadarajah, McGill University

Read More: https://arc.aiaa.org/doi/abs/10.2514/6.2017-3290
*/

/* SOURCE:
https://github.com/fanav/Verification_MMS_Tools
*/

extern double MS_1(int qnt_type, int qnt_id, double x, double y);
extern double MS_2(int qnt_type, int qnt_id, double x, double y);


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


// Main routine
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
    
//                            RANS-SA CASES                                    // 

    // MS-1: RANS-SA Verification case, original SA model
    case 1:  qnt = MS_1(qnt_type,qnt_id,x,y); break;        

    // MS-2: RANS-SA Verification case, modified SA model
    case 2:  qnt = MS_2(qnt_type,qnt_id,x,y); break;  
    
    // if nothing is chosen, exit the program.
    default:    printf("No MS is chosen! exit...\n"); exit(1);    
  } 
  
  return(qnt);
  
} 




//                            RANS-SA CASES                                    // 
double MS_1(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution forcing functions and quantities      
      
    // MS-1: RANS-SA Verification case, original SA model
      
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
  double dNTLdx,dNTLdy,d2NTLdxx,d2NTLdyy; //double d2NTLdxy; 
  double dPdx,dPdy,d2Pdxx,d2Pdyy;
  double Q_RHO_INV,Q_U_INV,Q_V_INV,Q_E_INV,Q_NTL_INV;
  double Q_RHO_VIS,Q_U_VIS,Q_V_VIS,Q_E_VIS,Q_NTL_VIS;  
  double Q_RHO,Q_U,Q_V,Q_E,Q_NTL;
  double Q_NTL_SRC_PRODUCT,Q_NTL_SRC_DESTRUCT,Q_NTL_SRC_DISTRIB,Q_NTL_SRC_CONSERV;   
  
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

          Q_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          Q_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          Q_RHO = Q_RHO_INV + Q_RHO_VIS; 
          qnt    = Q_RHO;
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
          
          Q_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = Q_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          Q_RHO_VIS = 0.;
          qnt    = Q_RHO_VIS;
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

          Q_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          Q_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          Q_U = Q_U_INV + Q_U_VIS; 
          qnt    = Q_U;
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

          Q_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U; 
          qnt    = Q_U_INV;
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
                 
          Q_U_VIS = ((2.0L/3.0L)*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(-(2*d2Udxx - d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) - (2*dUdx - dVdy)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(-(d2Udyy + d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) - (dUdy + dVdx)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(2*(2*dUdx - dVdy)*(dNTLdx*RHO + dRHOdx*NTL) + 3*(dUdy + dVdx)*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2))/pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2);           
          qnt    = Q_U_VIS;
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

          Q_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          Q_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          Q_V = Q_V_INV + Q_V_VIS; 
          qnt    = Q_V;
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
         
          Q_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V; 
          qnt    = Q_V_INV;
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
          
          Q_V_VIS = (-(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udxy + d2Vdxx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdy + dVdx)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (2.0L/3.0L)*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((d2Udxy - 2*d2Vdyy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (dUdx - 2*dVdy)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) + (mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(-2*(dUdx - 2*dVdy)*(dNTLdy*RHO + dRHOdy*NTL) + 3*(dUdy + dVdx)*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2))/pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2); 
          
          qnt    = Q_V_VIS;
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
          RHO    = MS_qnt_2D(MS_no,qt_Q,qi_Rho, x,y);
          U      = MS_qnt_2D(MS_no,qt_Q,qi_U,   x,y);
          V      = MS_qnt_2D(MS_no,qt_Q,qi_V,   x,y);
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

          Q_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          Q_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          Q_E = Q_E_INV + Q_E_VIS; 
          qnt    = Q_E;
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

          Q_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;
          qnt    = Q_E_INV;
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
          
          
          Q_E_VIS = ((2.0L/3.0L)*Pr*Prt*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(-dUdx*(2*dUdx - dVdy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + dVdy*(dUdx - 2*dVdy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) - (2*d2Udxx - d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*U + (d2Udxy - 2*d2Vdyy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*V + (dUdx - 2*dVdy)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)*V - (2*dUdx - dVdy)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)*U) - Pr*Prt*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(dUdy*(dUdy + dVdx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + dVdx*(dUdy + dVdx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4)) + (d2Udxy + d2Vdxx)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*V + (d2Udyy + d2Vdxy)*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*U + (dUdy + dVdx)*(4*dNTLdx*NTL*pow(RHO, 2) + 4*dRHOdx*pow(NTL, 2)*RHO + 3*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2)*V + (dUdy + dVdx)*(4*dNTLdy*NTL*pow(RHO, 2) + 4*dRHOdy*pow(NTL, 2)*RHO + 3*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)*U) + Pr*Prt*(mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)) + pow(NTL, 4)*pow(RHO, 4))*(-2*(dUdx - 2*dVdy)*(dNTLdy*RHO + dRHOdy*NTL)*V + 2*(2*dUdx - dVdy)*(dNTLdx*RHO + dRHOdx*NTL)*U + 3*(dUdy + dVdx)*(dNTLdx*RHO + dRHOdx*NTL)*V + 3*(dUdy + dVdx)*(dNTLdy*RHO + dRHOdy*NTL)*U)*pow(NTL, 2)*pow(RHO, 2) + Gamma*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*((Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*(-d2Edxx + 1.0*d2Udxx*U + 1.0*d2Vdxx*V + 1.0*pow(dUdx, 2) + 1.0*pow(dVdx, 2)) + (Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*(-d2Edyy + 1.0*d2Udyy*U + 1.0*d2Vdyy*V + 1.0*pow(dUdy, 2) + 1.0*pow(dVdy, 2)) + (-dEdx + 1.0*dUdx*U + 1.0*dVdx*V)*(4*Pr*dNTLdx*NTL*pow(RHO, 2) + 4*Pr*dRHOdx*pow(NTL, 2)*RHO + 3*Prt*mu*(dNTLdx*RHO + dRHOdx*NTL))*pow(NTL, 2)*pow(RHO, 2) + (-dEdy + 1.0*dUdy*U + 1.0*dVdy*V)*(4*Pr*dNTLdy*NTL*pow(RHO, 2) + 4*Pr*dRHOdy*pow(NTL, 2)*RHO + 3*Prt*mu*(dNTLdy*RHO + dRHOdy*NTL))*pow(NTL, 2)*pow(RHO, 2)) - 3*Gamma*((dNTLdx*RHO + dRHOdx*NTL)*(-dEdx + 1.0*dUdx*U + 1.0*dVdx*V) + (dNTLdy*RHO + dRHOdy*NTL)*(-dEdy + 1.0*dUdy*U + 1.0*dVdy*V))*(Pr*pow(NTL, 4)*pow(RHO, 4) + Prt*mu*(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3)))*pow(NTL, 2)*pow(RHO, 2))/(Pr*Prt*pow(pow(c_v1, 3)*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3), 2));
          qnt    = Q_E_VIS;
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

          Q_NTL_INV          = MS_qnt_2D(MS_no,qt_FFQ_inv,    qi_SA,x,y);
          Q_NTL_VIS          = MS_qnt_2D(MS_no,qt_FFQ_vis,    qi_SA,x,y);
          Q_NTL_SRC_PRODUCT  = MS_qnt_2D(MS_no,qt_FFQ_SA_prod,qi_SA,x,y);
          Q_NTL_SRC_DESTRUCT = MS_qnt_2D(MS_no,qt_FFQ_SA_dest,qi_SA,x,y);
          Q_NTL_SRC_DISTRIB  = MS_qnt_2D(MS_no,qt_FFQ_SA_dist,qi_SA,x,y);
          Q_NTL_SRC_CONSERV  = MS_qnt_2D(MS_no,qt_FFQ_SA_cons,qi_SA,x,y);
          
          Q_NTL = Q_NTL_INV + Q_NTL_VIS + Q_NTL_SRC_PRODUCT  + Q_NTL_SRC_DESTRUCT + Q_NTL_SRC_DISTRIB + Q_NTL_SRC_CONSERV; 
          qnt    = Q_NTL;
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

          Q_NTL_INV = dNTLdx*RHO*U + dNTLdy*RHO*V + dRHOdx*NTL*U + dRHOdy*NTL*V + dUdx*NTL*RHO + dVdy*NTL*RHO; 
          qnt    = Q_NTL_INV;
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
                    
          Q_NTL_VIS = -(d2NTLdxx*(mu + NTL*RHO) + d2NTLdyy*(mu + NTL*RHO) + dNTLdx*(dNTLdx*RHO + dRHOdx*NTL) + dNTLdy*(dNTLdy*RHO + dRHOdy*NTL))/sigma;          
          qnt = Q_NTL_VIS;
          break;
          
        // qnt_type=13: Production term in the SA source term
        case 13:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);          
          
          Q_NTL_SRC_PRODUCT = -c_b1*(-c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2)) + 1)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))*NTL*RHO;         
          qnt = Q_NTL_SRC_PRODUCT;
          break;       
          
        // qnt_type=14: Destruction term in the SA source term
        case 14:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);           
          
          Q_NTL_SRC_DESTRUCT = (-c_b1*c_t3*exp(-c_t4*pow(NTL, 2)*pow(RHO, 2)/pow(mu, 2))/pow(kappa, 2) + c_w1*pow((pow(c_w3, 6) + 1)/(pow(c_w3, 6) + pow(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6)), 1.0L/6.0L)*(c_w2*(pow(fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2))))), 6) - fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))) + fmin(r_lim, NTL/(pow(d, 2)*pow(kappa, 2)*(fabs(dUdy - dVdx) + (1 - NTL*RHO/(mu*(1 + pow(NTL, 4)*pow(RHO, 4)/(pow(mu, 4)*(pow(c_v1, 3) + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))))))*NTL/(pow(d, 2)*pow(kappa, 2)))))))*pow(NTL, 2)*RHO/pow(d, 2);           
          qnt = Q_NTL_SRC_DESTRUCT;
          break;  

        // qnt_type=15: Distribtuion term in the SA source term 
        case 15:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);
          
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                                 
          
          Q_NTL_SRC_DISTRIB = -c_b2*(pow(dNTLdx, 2) + pow(dNTLdy, 2))*RHO/sigma;           
          qnt = Q_NTL_SRC_DISTRIB;
          break;            
          
        // qnt_type=13: Conservation term in the SA source term 
        case 16:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
        
          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                       
                    
          Q_NTL_SRC_CONSERV = (dNTLdx*dRHOdx + dNTLdy*dRHOdy)*(mu/RHO + NTL)/sigma;
          qnt = Q_NTL_SRC_CONSERV;
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

double MS_2(int qnt_type, int qnt_id, double x, double y)
{
/*

      Manufactured solution forcing functions and quantities      
      
    // MS-2: RANS-SA Verification case, modified SA model
      
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
  double dNTLdx,dNTLdy,d2NTLdxx,d2NTLdyy; //double d2NTLdxy; 
  double dPdx,dPdy,d2Pdxx,d2Pdyy;
  double Q_RHO_INV,Q_U_INV,Q_V_INV,Q_E_INV,Q_NTL_INV;
  double Q_RHO_VIS,Q_U_VIS,Q_V_VIS,Q_E_VIS,Q_NTL_VIS;  
  double Q_RHO,Q_U,Q_V,Q_E,Q_NTL;  
  double Q_NTL_SRC_PRODUCT,Q_NTL_SRC_DESTRUCT,Q_NTL_SRC_DISTRIB,Q_NTL_SRC_CONSERV; 
  
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

          Q_RHO_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Cont,x,y);
          Q_RHO_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Cont,x,y);
          
          Q_RHO = Q_RHO_INV + Q_RHO_VIS; 
          qnt    = Q_RHO;
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
          
          Q_RHO_INV = dRHOdx*U + dRHOdy*V + dUdx*RHO + dVdy*RHO; 
          qnt    = Q_RHO_INV;
          break;          
          
        // qnt_type=12: viscous component of the forcing function           
        case 12:

          Q_RHO_VIS = 0.;
          qnt    = Q_RHO_VIS;
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

          Q_U_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momx,x,y);
          Q_U_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momx,x,y);
          
          Q_U = Q_U_INV + Q_U_VIS; 
          qnt    = Q_U;
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

          Q_U_INV = dPdx + dRHOdx*pow(U, 2) + dRHOdy*U*V + 2*dUdx*RHO*U + dUdy*RHO*V + dVdy*RHO*U; 
          qnt    = Q_U_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          d2Udxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_U,x,y);          
          d2Udyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_U,x,y);          
  
          d2Vdxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_V,x,y);         
                 
          Q_U_VIS = -1.0L/3.0L*mu*(4*d2Udxx + 3*d2Udyy + d2Vdxy);
          qnt    = Q_U_VIS;
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

          Q_V_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Momy,x,y);
          Q_V_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Momy,x,y);
          
          Q_V = Q_V_INV + Q_V_VIS; 
          qnt    = Q_V;
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
         
          Q_V_INV = dPdy + dRHOdx*U*V + dRHOdy*pow(V, 2) + dUdx*RHO*V + dVdx*RHO*U + 2*dVdy*RHO*V;
          qnt    = Q_V_INV;
          break;    
          
        // qnt_type=12: viscous component of the forcing function 
        case 12:
          d2Udxy = MS_qnt_2D(MS_no,qt_d2Qdxy,qi_U,x,y);          

          d2Vdxx = MS_qnt_2D(MS_no,qt_d2Qdxx,qi_V,x,y);                          
          d2Vdyy = MS_qnt_2D(MS_no,qt_d2Qdyy,qi_V,x,y);
                  
          Q_V_VIS = -1.0L/3.0L*mu*(d2Udxy + 3*d2Vdxx + 4*d2Vdyy);
          qnt    = Q_V_VIS;
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

          Q_E_INV = MS_qnt_2D(MS_no,qt_FFQ_inv,qi_Ener,x,y);
          Q_E_VIS = MS_qnt_2D(MS_no,qt_FFQ_vis,qi_Ener,x,y);
          
          Q_E = Q_E_INV + Q_E_VIS; 
          qnt    = Q_E;
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

          Q_E_INV = dUdx*(E*RHO + P) + dVdy*(E*RHO + P) + (dEdx*RHO + dPdx + dRHOdx*E)*U + (dEdy*RHO + dPdy + dRHOdy*E)*V;          
          qnt    = Q_E_INV;
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
          
          
          Q_E_VIS = (1.0L/3.0L)*mu*(Pr*(-2*dUdx*(2*dUdx - dVdy) - 3*dUdy*(dUdy + dVdx) - 3*dVdx*(dUdy + dVdx) + 2*dVdy*(dUdx - 2*dVdy) - 2*(2*d2Udxx - d2Vdxy)*U - 3*(d2Udxy + d2Vdxx)*V + 2*(d2Udxy - 2*d2Vdyy)*V - 3*(d2Udyy + d2Vdxy)*U) + 3*Gamma*(-d2Edxx + d2Udxx*U + d2Vdxx*V + pow(dUdx, 2) + pow(dVdx, 2)) + 3*Gamma*(-d2Edyy + d2Udyy*U + d2Vdyy*V + pow(dUdy, 2) + pow(dVdy, 2)))/Pr; 
          qnt    = Q_E_VIS;
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

          Q_NTL_INV          = MS_qnt_2D(MS_no,qt_FFQ_inv,    qi_SA,x,y);
          Q_NTL_VIS          = MS_qnt_2D(MS_no,qt_FFQ_vis,    qi_SA,x,y);
          Q_NTL_SRC_PRODUCT  = MS_qnt_2D(MS_no,qt_FFQ_SA_prod,qi_SA,x,y);
          Q_NTL_SRC_DESTRUCT = MS_qnt_2D(MS_no,qt_FFQ_SA_dest,qi_SA,x,y);
          Q_NTL_SRC_DISTRIB  = MS_qnt_2D(MS_no,qt_FFQ_SA_dist,qi_SA,x,y);
          Q_NTL_SRC_CONSERV  = MS_qnt_2D(MS_no,qt_FFQ_SA_cons,qi_SA,x,y);
          
          Q_NTL = Q_NTL_INV + Q_NTL_VIS + Q_NTL_SRC_PRODUCT  + Q_NTL_SRC_DESTRUCT + Q_NTL_SRC_DISTRIB + Q_NTL_SRC_CONSERV; 
          qnt    = Q_NTL;
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

          Q_NTL_INV = dNTLdx*RHO*U + dNTLdy*RHO*V + dRHOdx*NTL*U + dRHOdy*NTL*V + dUdx*NTL*RHO + dVdy*NTL*RHO; 
          qnt    = Q_NTL_INV;
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
          
          Q_NTL_VIS = -(dNTLdx*((c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3))*(dNTLdx*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*RHO + dRHOdx*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*NTL + 3*(dNTLdx*RHO + dRHOdx*NTL)*pow(NTL, 3)*pow(RHO, 3)) + 3*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(dNTLdx*RHO + dRHOdx*NTL)*pow(NTL, 3)*pow(RHO, 3)) + dNTLdy*((c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3))*(dNTLdy*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*RHO + dRHOdy*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*NTL + 3*(dNTLdy*RHO + dRHOdy*NTL)*pow(NTL, 3)*pow(RHO, 3)) + 3*(c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*(dNTLdy*RHO + dRHOdy*NTL)*pow(NTL, 3)*pow(RHO, 3)) + (d2NTLdxx + d2NTLdyy)*(c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3))*(mu*(c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3)) + (c_n1*pow(mu, 3) + pow(NTL, 3)*pow(RHO, 3))*NTL*RHO))/(sigma*pow(c_n1*pow(mu, 3) - pow(NTL, 3)*pow(RHO, 3), 2)); 
          qnt = Q_NTL_VIS;
          break;
          
        // qnt_type=13: Production term in the SA source term
        case 13:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dUdy     = MS_qnt_2D(MS_no,qt_dQdy,qi_U,x,y);
          
          dVdx     = MS_qnt_2D(MS_no,qt_dQdx,qi_V,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);          
          
          Q_NTL_SRC_PRODUCT = -c_b1*(-c_t3 + 1)*NTL*RHO*fabs(dUdy - dVdx);  
          qnt = Q_NTL_SRC_PRODUCT;
          break;       
          
        // qnt_type=14: Destruction term in the SA source term
        case 14:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          NTL      = MS_qnt_2D(MS_no,qt_Q,qi_Ntl,x,y);
          
          d        = MS_qnt_2D(MS_no,qt_Q,qi_Wall_d,x,y);            
          
          Q_NTL_SRC_DESTRUCT = -c_w1*pow(NTL, 2)*RHO/pow(d, 2); 
          qnt = Q_NTL_SRC_DESTRUCT;
          break;  

        // qnt_type=15: Distribtuion term in the SA source term 
        case 15:
          RHO      = MS_qnt_2D(MS_no,qt_Q,qi_Rho,x,y);

          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                       
          
          Q_NTL_SRC_DISTRIB = -c_b2*(pow(dNTLdx, 2) + pow(dNTLdy, 2))*RHO/sigma;
          qnt = Q_NTL_SRC_DISTRIB;
          break;            
          
        // qnt_type=13: Conservation term in the SA source term 
        case 16:
          RHO      = MS_qnt_2D(MS_no,qt_Q,   qi_Rho,x,y);
          dRHOdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Rho,x,y);
          dRHOdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Rho,x,y);
          
          NTL      = MS_qnt_2D(MS_no,qt_Q,   qi_Ntl,x,y);
          dNTLdx   = MS_qnt_2D(MS_no,qt_dQdx,qi_Ntl,x,y);
          dNTLdy   = MS_qnt_2D(MS_no,qt_dQdy,qi_Ntl,x,y);                       
          
          
          Q_NTL_SRC_CONSERV = (dNTLdx*dRHOdx + dNTLdy*dRHOdy)*(mu/RHO + (c_n1 + pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3))*NTL/(c_n1 - pow(NTL, 3)*pow(RHO, 3)/pow(mu, 3)))/sigma; 
          qnt = Q_NTL_SRC_CONSERV;
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


