#include "RadiativeCooling.h"



m_radcool = new RadiativeCooling(outpath);

float Tcmb = 2.725f;//just in case someone set the CMB temp to zero
m_radcool->setTCMB(Tcmb);//Set CMB Temp
m_radcool->readCloudyScaleFactor((RAD_T)aa);// find rough a value

    // sph volume estimate
    POSVEL_T Vi = mass[nn]/rho[nn];



    RAD_T rhoi = rho[nn]*m_hubble*m_hubble*omega_cb*(1.0f/(aa*aa*aa));
    RAD_T mui  = mu[nn];
    RAD_T Ti   = CONV_U_TO_T*UU_CGS*scale_uu*mui*uu[nn]*aa*aa;
    RAD_T Zi   = zmet[nn];
    RAD_T Yi   = yhe[nn];
    RAD_T Xi = (1.0-Yi-Zi);
    RAD_T nHi  = rhoi*RHO_CGS*Xi*INV_MH;
    RAD_T nHIi = 0.0;
    RAD_T nei = 0.0;
#ifdef UV_BACKGROUND
    int iter = 0;
#ifndef UVCONVERT
    rhoi *= RHO_CGS;
#endif



RAD_T lambda = (*m_radcool)(Ti, rhoi, Zi, Yi, (RAD_T)aa, mui, &iter, false, &nHIi, &nei);
nei *= nHi

// the else is
// - check this is consistent 
    nei = (1.0f - mui*(1.0f-0.75f*PRIMORDIAL_Y))/((1.0f-PRIMORDIAL_Y)*mui);//formula for ne/nH
    nei *= nHi; // ne/nH -> ne


