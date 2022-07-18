#include "../eos/eos.hpp"
#include "../index_def.hpp"
#include "hllc.hpp"
#include <cmath>
#include <algorithm>

void hllc( double *valsL,
          double *valsR,
          double *fluxs,
          int IMP,
          double &gamma){
    /**
        rhoL, rhoR: density
        pL, pR:     pressure
        uL, uR:     speed
        valsL:      array of double storing the values on left. Quantities are <rho, rho*v1, rho*v2, rho*v3, E> (may not in this order)
        gamma:      hydro gamma
        INDRHO:     the ordering of variables in valsL and valsR (density)
        INDV1:      the ordering of variables in valsL and valsR (v1)
        INDV2:      the ordering of variables in valsL and valsR (v2)
        INDV3:      the ordering of variables in valsL and valsR (v3)
        INDE :      the ordering of variables in valsL and valsR (E)
    **/
    // step 1: wave speed estimates, Toro 10.6
    double rhoL = valsL[IDN];
    double rhoR = valsR[IDN];
    double EL = valsL[IEN];
    double ER = valsR[IEN];
    double pL = pres(valsL[IDN], valsL[IEN], valsL[IM1], valsL[IM2], valsL[IM3], gamma);
    double pR = pres(valsR[IDN], valsR[IEN], valsR[IM1], valsR[IM2], valsR[IM3], gamma);
    double uL = valsL[IMP] / rhoL;
    double uR = valsR[IMP] / rhoR;
    double aL = soundspeed(rhoL, pL, gamma);
    double aR = soundspeed(rhoR, pR, gamma);
    double rhobar = 0.5 * (rhoL + rhoR);
    double abar   = 0.5 * (aL + aR);
    double ppvrs  = 0.5 * (pL + pR) - 0.5 * (uR - uL) * rhobar * abar;
    double pstar  = std::max(0., ppvrs);
    // step 2: wave speed estimates
    double qL, qR;
    if (pstar <= pL){
        qL = 1;
    }
    else{
        qL = sqrt(1 + (gamma + 1) / (2 * gamma) * (pstar / pL - 1));
    }

    if (pstar <= pR){
        qR = 1;
    }
    else{
        qR = sqrt(1 + (gamma + 1) / (2 * gamma) * (pstar / pR - 1));
    }

    double sL = uL - aL * qL;
    double sR = uR + aR * qR;

    double sstar = (pR - pL + rhoL * uL * (sL - uL) - rhoR * uR * (sR - uR)) / (rhoL * (sL - uL) - rhoR * (sR - uR));     // 10.70
    // step 3: HLLC flux
    for (int val_ind = 0; val_ind < 5; val_ind ++){
        double flux_L = valsL[val_ind] * uL;
        double flux_R = valsR[val_ind] * uR;
        if (val_ind == IMP){
            flux_L += pL;
            flux_R += pR;
        }
        if (val_ind == IEN){
            flux_L += pL * uL;
            flux_R += pR * uR;
        }
        //cout << "flux: " << flux_L << '\t' << flux_R << endl;
        if      (0 <= sL){
            fluxs[val_ind] = flux_L;
        }
        else if (sL <= 0 <= sstar){
            double UsKM;
            if (val_ind == IDN)     { UsKM = 1; }
            else if (val_ind == IEN){ UsKM = EL / rhoL + (sstar - uL) * (sstar + pL / (rhoL * (sL - uL))); }
            else if (val_ind == IMP){ UsKM = sstar; }
            else                    { UsKM = valsL[val_ind] / valsL[IDN]; }
            double ustarL = rhoL * (sL - uL) / (sL - sstar);
            fluxs[val_ind] = flux_L + sL * (ustarL - uL);
        }
        else if (sstar <= 0 <= sR){
            double UsKM;
            if (val_ind == IDN)     { UsKM = 1; }
            else if (val_ind == IEN){ UsKM = ER / rhoR + (sstar - uR) * (sstar + pR / (rhoR * (sR - uR))); }
            else if (val_ind == IMP){ UsKM = sstar; }
            else                    { UsKM = valsR[val_ind] / valsR[IDN]; }
            double ustarR = rhoR * (sR - uR) / (sR - sstar);
            fluxs[val_ind] = flux_R + sR * (ustarR - uR);
        }
        else {     //  if (0 >= sR )
            fluxs[val_ind] = flux_R;
        }
    }
}
