#include "../eos/eos.hpp"
#include "../index_def.hpp"
#include <cmath>


void hll_dust( double *valsL,
          double *valsR,
          double *fluxs,
          int IMP,
          double &gamma){

    // calculate pressure and velocities
    double uL = valsL[IMP] / valsL[IDN];
    double uR = valsR[IMP] / valsR[IDN];
    // step 1: calculate flux
    for (int val_ind = 0; val_ind < 4; val_ind ++){
        double flux_L = valsL[val_ind] * uL;
        double flux_R = valsR[val_ind] * uR;

        if      (0 <= uL){ fluxs[val_ind] = flux_L; }
        else if (0 >= uR ){ fluxs[val_ind] = flux_R; }
        else{
            fluxs[val_ind] = (uR * flux_L - uL * flux_R + uL * uR * (valsR[val_ind] - valsL[val_ind])) / (uR - uL);
        }
    }
}

