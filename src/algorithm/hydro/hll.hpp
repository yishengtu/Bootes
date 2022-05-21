#ifndef HLL_HPP_
#define HLL_HPP_

#include "../eos/eos.hpp"
#include <cmath>


void hll( BootesArray<float> &valsL,
          BootesArray<float> &valsR,
          BootesArray<float> &fluxs,
          int IMP,
          float &gamma){

    // calculate pressure and velocities
    float pL = pres(valsL(IDN), valsL(IEN), valsL(IM1), valsL(IM2), valsL(IM3), gamma);
    float pR = pres(valsR(IDN), valsR(IEN), valsR(IM1), valsR(IM2), valsR(IM3), gamma);
    float uL = valsL(IMP) / valsL(IDN);
    float uR = valsR(IMP) / valsR(IDN);
    // step 1: wave speed estimates, Toro 10.5.1
    float aL = soundspeed(valsL(IDN), pL, gamma);
    float aR = soundspeed(valsR(IDN), pR, gamma);
    float sL = uL - aL;
    float sR = uR + aR;
    // step 2: calculate flux
    for (int val_ind = 0; val_ind < fluxs.shape()[0]; val_ind ++){
        float flux_L = valsL(val_ind) * uL;
        float flux_R = valsR(val_ind) * uR;
        if (val_ind == IMP){
            flux_L += pL;
            flux_R += pR;
        }
        if (val_ind == IEN){
            flux_L += pL * uL;
            flux_R += pR * uR;
        }
        //cout << "flux: " << flux_L << '\t' << flux_R << endl;
        if      (0 <= sL){ fluxs(val_ind) = flux_L; }
        else if (0 >= sR ){ fluxs(val_ind) = flux_R; }
        else{
            fluxs(val_ind) = (sR * flux_L - sL * flux_R + sL * sR * (valsR(val_ind) - valsL(val_ind))) / (sR - sL);
        }
    }
}


#endif // HLL_HPP_
