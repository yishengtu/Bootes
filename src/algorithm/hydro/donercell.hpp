#ifndef DONERCELL_HPP_
#define DONERCELL_HPP_

#include "../eos/eos.hpp"
#include <cmath>


void donercell(double &rhoL, double &pL, double &uL,
          double &rhoR, double &pR, double &uR,
          BootesArray<double> &valsL,
          BootesArray<double> &valsR,
          BootesArray<double> &fluxs,
          int MOMPRESIND,
          double &gamma){

    // step 1: wave speed estimates, Toro 10.5.1
    double aL = soundspeed(rhoL, pL, gamma);
    double aR = soundspeed(rhoR, pR, gamma);
    double sL = uL - aL;
    double sR = uR + aR;
    // step 2: calculate flux
    for (int val_ind = 0; val_ind < fluxs.shape()[0]; val_ind ++){
        double flux_L = valsL(val_ind) * uL;
        double flux_R = valsR(val_ind) * uR;
        if (val_ind == MOMPRESIND){
            flux_L += pL;
            flux_R += pR;
        }
        fluxs(val_ind) = flux_L + flux_R;
    }
}


#endif // HLL_HPP_
