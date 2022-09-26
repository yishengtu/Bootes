#include "../BootesArray.hpp"
#include "../index_def.hpp"


void sph_polar_pole_boundary_condition_x2i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                                int &x2s, int &x2l, int &ng2,
                                                                int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //#pragma omp parallel for collapse(3)
    #pragma acc parallel loop collapse(3) default (present) firstprivate(x1s, x1l, ng1, x2s, x2l, ng2, x3s, x3l, ng3)
    for (int gind2 = 0; gind2 < ng2; gind2 ++){
        for (int kk = x3s; kk < x3l; kk++){
            for (int ii = x1s; ii < x1l; ii++){
                quan(IDN, kk, x2s - 1 - gind2, ii) = quan(IDN, kk, x2s + gind2, ii);
                quan(IM1, kk, x2s - 1 - gind2, ii) = quan(IM1, kk, x2s + gind2, ii);
                quan(IM2, kk, x2s - 1 - gind2, ii) = - quan(IM2, kk, x2s + gind2, ii);
                quan(IM3, kk, x2s - 1 - gind2, ii) = - quan(IM3, kk, x2s + gind2, ii);
                quan(IEN, kk, x2s - 1 - gind2, ii) = quan(IEN, kk, x2s + gind2, ii);
            }
        }
    }
}


void sph_polar_pole_boundary_condition_x2o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                                int &x2s, int &x2l, int &ng2,
                                                                int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //#pragma omp parallel for collapse(3)
    #pragma acc parallel loop collapse(3) default (present) firstprivate(x1s, x1l, ng1, x2s, x2l, ng2, x3s, x3l, ng3)
    for (int gind2 = 0; gind2 < ng2; gind2 ++){
        for (int kk = x3s; kk < x3l; kk++){
            for (int ii = x1s; ii < x1l; ii++){
                quan(IEN, kk, x2l + gind2, ii)     = quan(IEN, kk, x2l - (gind2 + 1), ii);
                quan(IM1, kk, x2l + gind2, ii)     = quan(IM1, kk, x2l - (gind2 + 1), ii);
                quan(IM2, kk, x2l + gind2, ii)     = - quan(IM2, kk, x2l - (gind2 + 1), ii);
                quan(IM3, kk, x2l + gind2, ii)     = - quan(IM3, kk, x2l - (gind2 + 1), ii);
                quan(IDN, kk, x2l + gind2, ii)     = quan(IDN, kk, x2l - (gind2 + 1), ii);
            }
        }
    }
}

