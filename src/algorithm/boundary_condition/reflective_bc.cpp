#include "../BootesArray.hpp"
#include "../index_def.hpp"
#include "reflective_bc.hpp"


void reflective_boundary_condition_x1i(BootesArray<double> &quan, int x1s, int x1l, int ng1,
                                                                int x2s, int x2l, int ng2,
                                                                int x3s, int x3l, int ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //#pragma omp parallel for collapse(3) schedule (static)
    #pragma acc parallel loop collapse(3) default (present)
    for (int gind1 = 0; gind1 < ng1; gind1 ++){
        for (int kk = x3s; kk < x3l; kk++){
            for (int jj = x2s; jj < x2l; jj++){
                quan(IDN, kk, jj, x1s - 1 - gind1) = quan(IDN, kk, jj, x1s + gind1);
                quan(IM1, kk, jj, x1s - 1 - gind1) = - quan(IM1, kk, jj, x1s + gind1);
                quan(IM2, kk, jj, x1s - 1 - gind1) = quan(IM2, kk, jj, x1s + gind1);
                quan(IM3, kk, jj, x1s - 1 - gind1) = quan(IM3, kk, jj, x1s + gind1);
                quan(IEN, kk, jj, x1s - 1 - gind1) = quan(IEN, kk, jj, x1s + gind1);
            }
        }
    }
}


void reflective_boundary_condition_x1o(BootesArray<double> &quan, int x1s, int x1l, int ng1,
                                                                int x2s, int x2l, int ng2,
                                                                int x3s, int x3l, int ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //#pragma omp parallel for collapse(3) schedule (static)
    #pragma acc parallel loop collapse(3) default (present)
    for (int gind1 = 0; gind1 < ng1; gind1 ++){
        for (int kk = x3s; kk < x3l; kk++){
            for (int jj = x2s; jj < x2l; jj++){
                quan(IDN, kk, jj, x1l + gind1)     = quan(IDN, kk, jj, x1l - (gind1 + 1));
                quan(IM1, kk, jj, x1l + gind1)     = - quan(IM1, kk, jj, x1l - (gind1 + 1));
                quan(IM2, kk, jj, x1l + gind1)     = quan(IM2, kk, jj, x1l - (gind1 + 1));
                quan(IM3, kk, jj, x1l + gind1)     = quan(IM3, kk, jj, x1l - (gind1 + 1));
                quan(IEN, kk, jj, x1l + gind1)     = quan(IEN, kk, jj, x1l - (gind1 + 1));
            }
        }
    }
}


void reflective_boundary_condition_x2i(BootesArray<double> &quan, int x1s, int x1l, int ng1,
                                                                int x2s, int x2l, int ng2,
                                                                int x3s, int x3l, int ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //#pragma omp parallel for collapse(3) schedule (static)
    #pragma acc parallel loop collapse(3) default (present)
    for (int gind2 = 0; gind2 < ng2; gind2 ++){
        for (int kk = x3s; kk < x3l; kk++){
            for (int ii = x1s; ii < x1l; ii++){
                quan(IDN, kk, x2s - 1 - gind2, ii) = quan(IDN, kk, x2s + gind2, ii);
                quan(IM1, kk, x2s - 1 - gind2, ii) = quan(IM1, kk, x2s + gind2, ii);
                quan(IM2, kk, x2s - 1 - gind2, ii) = - quan(IM2, kk, x2s + gind2, ii);
                quan(IM3, kk, x2s - 1 - gind2, ii) = quan(IM3, kk, x2s + gind2, ii);
                quan(IEN, kk, x2s - 1 - gind2, ii) = quan(IEN, kk, x2s + gind2, ii);
            }
        }
    }
}


void reflective_boundary_condition_x2o(BootesArray<double> &quan, int x1s, int x1l, int ng1,
                                                                  int x2s, int x2l, int ng2,
                                                                  int x3s, int x3l, int ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //#pragma omp parallel for collapse(3) schedule (static)
    #pragma acc parallel loop collapse(3) default (present)
    for (int gind2 = 0; gind2 < ng2; gind2 ++){
        for (int kk = x3s; kk < x3l; kk++){
            for (int ii = x1s; ii < x1l; ii++){
                quan(IDN, kk, x2l + gind2, ii)     = quan(IDN, kk, x2l - (gind2 + 1), ii);
                quan(IM1, kk, x2l + gind2, ii)     = quan(IM1, kk, x2l - (gind2 + 1), ii);
                quan(IM2, kk, x2l + gind2, ii)     = - quan(IM2, kk, x2l - (gind2 + 1), ii);
                quan(IM3, kk, x2l + gind2, ii)     = quan(IM3, kk, x2l - (gind2 + 1), ii);
                quan(IEN, kk, x2l + gind2, ii)     = quan(IEN, kk, x2l - (gind2 + 1), ii);
            }
        }
    }
}


void reflective_boundary_condition_x3i(BootesArray<double> &quan, int x1s, int x1l, int ng1,
                                                                  int x2s, int x2l, int ng2,
                                                                  int x3s, int x3l, int ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //#pragma omp parallel for collapse(3) schedule (static)
    #pragma acc parallel loop collapse(3) default (present)
    for (int gind3 = 0; gind3 < ng3; gind3 ++){
        for (int jj = x2s; jj < x2l; jj++){
            for (int ii = x1s; ii < x1l; ii++){
                quan(IDN, x3s - 1 - gind3, jj, ii) = quan(IDN, x3s + gind3, jj, ii);
                quan(IM1, x3s - 1 - gind3, jj, ii) = quan(IM1, x3s + gind3, jj, ii);
                quan(IM2, x3s - 1 - gind3, jj, ii) = quan(IM2, x3s + gind3, jj, ii);
                quan(IM3, x3s - 1 - gind3, jj, ii) = - quan(IM3, x3s + gind3, jj, ii);
                quan(IEN, x3s - 1 - gind3, jj, ii) = quan(IEN, x3s + gind3, jj, ii);
            }
        }
    }
}


void reflective_boundary_condition_x3o(BootesArray<double> &quan, int x1s, int x1l, int ng1,
                                                                  int x2s, int x2l, int ng2,
                                                                  int x3s, int x3l, int ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //#pragma omp parallel for collapse(3) schedule (static)
    #pragma acc parallel loop collapse(3) default (present)
    for (int gind3 = 0; gind3 < ng3; gind3 ++){
        for (int jj = x2s; jj < x2l; jj++){
            for (int ii = x1s; ii < x1l; ii++){
                quan(IDN, x3l + gind3, jj, ii)     = quan(IDN, x3l - (gind3 + 1), jj, ii);
                quan(IM1, x3l + gind3, jj, ii)     = quan(IM1, x3l - (gind3 + 1), jj, ii);
                quan(IM2, x3l + gind3, jj, ii)     = quan(IM2, x3l - (gind3 + 1), jj, ii);
                quan(IM3, x3l + gind3, jj, ii)     = - quan(IM3, x3l - (gind3 + 1), jj, ii);
                quan(IEN, x3l + gind3, jj, ii)     = quan(IEN, x3l - (gind3 + 1), jj, ii);
            }
        }
    }
}


