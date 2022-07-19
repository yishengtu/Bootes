#include "../eos/eos.hpp"
#include "../index_def.hpp"
#include "hlle.hpp"
#include <cmath>
#include <algorithm>

void hlle(double *valsL,
          double *valsR,
          double *fluxs,
          int IMP,
          double &gamma){

    // calculate pressure and velocities
    double pL = pres(valsL[IDN], valsL[IEN], valsL[IM1], valsL[IM2], valsL[IM3], gamma);
    double pR = pres(valsR[IDN], valsR[IEN], valsR[IM1], valsR[IM2], valsR[IM3], gamma);
    double vL = valsL[IMP] / valsL[IDN];
    double vR = valsR[IMP] / valsR[IDN];
    // step 1: wave speed estimates
    double cL = soundspeed(valsL[IDN], pL, gamma);
    double cR = soundspeed(valsR[IDN], pR, gamma);
    double aL = std::min(vL - cL, vR - cR);
    double aR = std::max(vL + cL, vR + cR);
    double bp = std::max(aR, (double) 0);
    double bm = std::min(aL, (double) 0);
    double vxL = vL - bm;
    double vxR = vR - bp;
    // step 2: calculate flux
    for (int val_ind = 0; val_ind < 5; val_ind ++){
        double flux_L = valsL[val_ind] * vxL;
        double flux_R = valsR[val_ind] * vxR;
        if (val_ind == IMP){
            flux_L += pL;
            flux_R += pR;
        }
        if (val_ind == IEN){
            flux_L += pL * vL;
            flux_R += pR * vR;
        }
        double tmp = 0;
        if (bp != bm){
             tmp = 0.5 * (bm + bp) / (bp - bm);
        }
        fluxs[val_ind] = 0.5 * (flux_L + flux_R) + (flux_L - flux_R) * tmp;
    }
}


