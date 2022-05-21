#ifndef DONERCELL_HPP_
#define DONERCELL_HPP_

#include "../eos/eos.hpp"
#include <cmath>


void donercell(float &rhoL, float &pL, float &uL,
          float &rhoR, float &pR, float &uR,
          BootesArray<float> &valsL,
          BootesArray<float> &valsR,
          BootesArray<float> &fluxs,
          int MOMPRESIND,
          float &gamma){

    // step 1: wave speed estimates, Toro 10.5.1
    float aL = soundspeed(rhoL, pL, gamma);
    float aR = soundspeed(rhoR, pR, gamma);
    float sL = uL - aL;
    float sR = uR + aR;
    // step 2: calculate flux
    for (int val_ind = 0; val_ind < fluxs.shape()[0]; val_ind ++){
        float flux_L = valsL(val_ind) * uL;
        float flux_R = valsR(val_ind) * uR;
        if (val_ind == MOMPRESIND){
            flux_L += pL;
            flux_R += pR;
        }
        fluxs(val_ind) = flux_L + flux_R;
    }
}


#endif // HLL_HPP_
