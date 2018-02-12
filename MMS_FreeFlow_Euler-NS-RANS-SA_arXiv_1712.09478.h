#ifdef __cplusplus
extern “C”
#else
extern
#endif

extern double  MS_qnt_2D(int MS_no, int qnt_type, int qnt_id, double x, double y);

// In order to reduce the risk of bug, use enum types instead of integers
// to set the values of 'qnt_type' and 'qnt_id' arguments in calls to MS_qnt_2D.
// To deactivate enum definitions, use #define NO_ENUM_TO_CALL_MS_qnt_2D at a higher inclusion level.
#ifndef NO_ENUM_TO_CALL_MS_qnt_2D
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
#endif


