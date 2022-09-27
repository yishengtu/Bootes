#include "../../BootesArray.hpp"
#include "../../index_def.hpp"


void dust_reflective_boundary_condition_x1i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                                int &x2s, int &x2l, int &ng2,
                                                                int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    // #pragma omp parallel for collapse(4) schedule (static)
    int shape0 = quan.shape()[0];
    #pragma acc parallel loop collapse(4) default (present) firstprivate(x1s, x1l, ng1, x2s, x2l, ng2, x3s, x3l, ng3)
    for (int specIND = 0; specIND < shape0; specIND ++){
        for (int gind1 = 0; gind1 < ng1; gind1 ++){
            for (int kk = x3s; kk < x3l; kk++){
                for (int jj = x2s; jj < x2l; jj++){
                    quan(specIND, IDN, kk, jj, x1s - 1 - gind1) = quan(specIND, IDN, kk, jj, x1s + gind1);
                    quan(specIND, IM1, kk, jj, x1s - 1 - gind1) = - quan(specIND, IM1, kk, jj, x1s + gind1);
                    quan(specIND, IM2, kk, jj, x1s - 1 - gind1) = quan(specIND, IM2, kk, jj, x1s + gind1);
                    quan(specIND, IM3, kk, jj, x1s - 1 - gind1) = quan(specIND, IM3, kk, jj, x1s + gind1);
                }
            }
        }
    }
}


void dust_reflective_boundary_condition_x1o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                                int &x2s, int &x2l, int &ng2,
                                                                int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    // #pragma omp parallel for collapse(4) schedule (static)
    int shape0 = quan.shape()[0];
    #pragma acc parallel loop collapse(4) default (present) firstprivate(x1s, x1l, ng1, x2s, x2l, ng2, x3s, x3l, ng3)
    for (int specIND = 0; specIND < shape0; specIND ++){
        for (int gind1 = 0; gind1 < ng1; gind1 ++){
            for (int kk = x3s; kk < x3l; kk++){
                for (int jj = x2s; jj < x2l; jj++){
                    quan(specIND, IDN, kk, jj, x1l + gind1)     = quan(specIND, IDN, kk, jj, x1l - (gind1 + 1));
                    quan(specIND, IM1, kk, jj, x1l + gind1)     = - quan(specIND, IM1, kk, jj, x1l - (gind1 + 1));
                    quan(specIND, IM2, kk, jj, x1l + gind1)     = quan(specIND, IM2, kk, jj, x1l - (gind1 + 1));
                    quan(specIND, IM3, kk, jj, x1l + gind1)     = quan(specIND, IM3, kk, jj, x1l - (gind1 + 1));
                }
            }
        }
    }
}


void dust_reflective_boundary_condition_x2i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                                int &x2s, int &x2l, int &ng2,
                                                                int &x3s, int &x3l, int &ng3){
    int shape0 = quan.shape()[0];
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    // #pragma omp parallel for collapse(4) schedule (static)
    #pragma acc parallel loop collapse(4) default (present) firstprivate(x1s, x1l, ng1, x2s, x2l, ng2, x3s, x3l, ng3)
    for (int specIND = 0; specIND < shape0; specIND ++){
        for (int gind2 = 0; gind2 < ng2; gind2 ++){
            for (int kk = x3s; kk < x3l; kk++){
                for (int ii = x1s; ii < x1l; ii++){
                    quan(specIND, IDN, kk, x2s - 1 - gind2, ii) = quan(specIND, IDN, kk, x2s + gind2, ii);
                    quan(specIND, IM1, kk, x2s - 1 - gind2, ii) = quan(specIND, IM1, kk, x2s + gind2, ii);
                    quan(specIND, IM2, kk, x2s - 1 - gind2, ii) = - quan(specIND, IM2, kk, x2s + gind2, ii);
                    quan(specIND, IM3, kk, x2s - 1 - gind2, ii) = quan(specIND, IM3, kk, x2s + gind2, ii);
                }
            }
        }
    }
}


void dust_reflective_boundary_condition_x2o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                                int &x2s, int &x2l, int &ng2,
                                                                int &x3s, int &x3l, int &ng3){
    int shape0 = quan.shape()[0];
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    // #pragma omp parallel for collapse(4) schedule (static)
    #pragma acc parallel loop collapse(4) default (present) firstprivate(x1s, x1l, ng1, x2s, x2l, ng2, x3s, x3l, ng3)
    for (int specIND = 0; specIND < shape0; specIND ++){
        for (int gind2 = 0; gind2 < ng2; gind2 ++){
            for (int kk = x3s; kk < x3l; kk++){
                for (int ii = x1s; ii < x1l; ii++){
                    quan(specIND, IDN, kk, x2l + gind2, ii)     = quan(specIND, IDN, kk, x2l - (gind2 + 1), ii);
                    quan(specIND, IM1, kk, x2l + gind2, ii)     = quan(specIND, IM1, kk, x2l - (gind2 + 1), ii);
                    quan(specIND, IM2, kk, x2l + gind2, ii)     = - quan(specIND, IM2, kk, x2l - (gind2 + 1), ii);
                    quan(specIND, IM3, kk, x2l + gind2, ii)     = quan(specIND, IM3, kk, x2l - (gind2 + 1), ii);
                }
            }
        }
    }
}


void dust_reflective_boundary_condition_x3i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                                int &x2s, int &x2l, int &ng2,
                                                                int &x3s, int &x3l, int &ng3){
    int shape0 = quan.shape()[0];
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    // #pragma omp parallel for collapse(4) schedule (static)
    #pragma acc parallel loop collapse(4) default (present) firstprivate(x1s, x1l, ng1, x2s, x2l, ng2, x3s, x3l, ng3)
    for (int specIND = 0; specIND < shape0; specIND ++){
        for (int gind3 = 0; gind3 < ng3; gind3 ++){
            for (int jj = x2s; jj < x2l; jj++){
                for (int ii = x1s; ii < x1l; ii++){
                    quan(specIND, IDN, x3s - 1 - gind3, jj, ii) = quan(specIND, IDN, x3s + gind3, jj, ii);
                    quan(specIND, IM1, x3s - 1 - gind3, jj, ii) = quan(specIND, IM1, x3s + gind3, jj, ii);
                    quan(specIND, IM2, x3s - 1 - gind3, jj, ii) = quan(specIND, IM2, x3s + gind3, jj, ii);
                    quan(specIND, IM3, x3s - 1 - gind3, jj, ii) = - quan(specIND, IM3, x3s + gind3, jj, ii);
                }
            }
        }
    }
}


void dust_reflective_boundary_condition_x3o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                                int &x2s, int &x2l, int &ng2,
                                                                int &x3s, int &x3l, int &ng3){
    int shape0 = quan.shape()[0];
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    //  #pragma omp parallel for collapse(4) schedule (static)
    #pragma acc parallel loop collapse(4) default (present) firstprivate(x1s, x1l, ng1, x2s, x2l, ng2, x3s, x3l, ng3)
    for (int specIND = 0; specIND < shape0; specIND ++){
        for (int gind3 = 0; gind3 < ng3; gind3 ++){
            for (int jj = x2s; jj < x2l; jj++){
                for (int ii = x1s; ii < x1l; ii++){
                    quan(specIND, IDN, x3l + gind3, jj, ii) = quan(specIND, IDN, x3l - (gind3 + 1), jj, ii);
                    quan(specIND, IM1, x3l + gind3, jj, ii) = quan(specIND, IM1, x3l - (gind3 + 1), jj, ii);
                    quan(specIND, IM2, x3l + gind3, jj, ii) = quan(specIND, IM2, x3l - (gind3 + 1), jj, ii);
                    quan(specIND, IM3, x3l + gind3, jj, ii) = - quan(specIND, IM3, x3l - (gind3 + 1), jj, ii);
                }
            }
        }
    }
}


