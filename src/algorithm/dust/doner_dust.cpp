#include "../eos/eos.hpp"
#include "../index_def.hpp"
#include "doner_dust.hpp"
#include <cmath>


void doner_cell_dust( double *valsL,
          double *valsR,
          double *fluxs,
          int IMP,
          double &gamma){

    // calculate pressure and velocities
    double uL = valsL[IMP] / valsL[IDN];
    double uR = valsR[IMP] / valsR[IDN];
    // step 1: wave speed estimates, Toro 10.5.1, when no pressure, no sound speed. so s = u
    // double sL = uL;
    // double sR = uR;
    // step 2: calculate flux
    for (int val_ind = 0; val_ind < 4; val_ind ++){
        double flux_L = valsL[val_ind] * uL;
        double flux_R = valsR[val_ind] * uR;
        if      (uL >= 0 && uR >= 0){ fluxs[val_ind] = flux_L; }
        else if (uR <= 0 && uL <= 0){ fluxs[val_ind] = flux_R; }
        else if (uR >= 0 && uL <= 0){ fluxs[val_ind] = 0; }
        else{
            fluxs[val_ind] = flux_L + flux_R;
        }
    }
}

