#ifndef HLL_HPP_
#define HLL_HPP_

#include "../eos/eos.hpp"
#include <cmath>


void hll( BootesArray<double> &valsL,
          BootesArray<double> &valsR,
          BootesArray<double> &fluxs,
          int IMP,
          double &gamma){

    // calculate pressure and velocities
    double pL = pres(valsL(IDN), valsL(IEN), valsL(IM1), valsL(IM2), valsL(IM3), gamma);
    double pR = pres(valsR(IDN), valsR(IEN), valsR(IM1), valsR(IM2), valsR(IM3), gamma);
    double uL = valsL(IMP) / valsL(IDN);
    double uR = valsR(IMP) / valsR(IDN);
    // step 1: wave speed estimates, Toro 10.5.1
    double aL = soundspeed(valsL(IDN), pL, gamma);
    double aR = soundspeed(valsR(IDN), pR, gamma);
    double sL = uL - aL;
    double sR = uR + aR;
    // step 2: calculate flux
    for (int val_ind = 0; val_ind < fluxs.shape()[0]; val_ind ++){
        double flux_L = valsL(val_ind) * uL;
        double flux_R = valsR(val_ind) * uR;
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
