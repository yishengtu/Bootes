#include "../../BootesArray.hpp"
#include "../../index_def.hpp"
#include "standard_bc_dust.hpp"


void dust_standard_boundary_condition_x1i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    #pragma omp parallel for collapse(4)
    for (int specIND = 0; specIND < quan.shape()[0]; specIND ++){
        for (int valIND = 0; valIND < quan.shape()[1]; valIND ++){
            for (int gind1 = 0; gind1 < ng1; gind1 ++){
                for (int kk = x3s; kk < x3l; kk++){
                    for (int jj = x2s; jj < x2l; jj++){
                        quan(specIND, valIND, kk, jj, x1s - 1 - gind1) = quan(specIND, valIND, kk, jj, x1s + gind1);
                    }
                }
            }
        }
    }
}


void dust_standard_boundary_condition_x1o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    #pragma omp parallel for collapse(4)
    for (int specIND = 0; specIND < quan.shape()[0]; specIND ++){
        for (int valIND = 0; valIND < quan.shape()[1]; valIND ++){
            for (int gind1 = 0; gind1 < ng1; gind1 ++){
                for (int kk = x3s; kk < x3l; kk++){
                    for (int jj = x2s; jj < x2l; jj++){
                        quan(specIND, valIND, kk, jj, x1l + gind1)     = quan(specIND, valIND, kk, jj, x1l - (gind1 + 1));
                    }
                }
            }
        }
    }
}


void dust_standard_boundary_condition_x2i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    #pragma omp parallel for collapse(4)
    for (int specIND = 0; specIND < quan.shape()[0]; specIND ++){
        for (int valIND = 0; valIND < quan.shape()[1]; valIND ++){
            for (int gind2 = 0; gind2 < ng2; gind2 ++){
                for (int kk = x3s; kk < x3l; kk++){
                    for (int ii = x1s; ii < x1l; ii++){
                        quan(specIND, valIND, kk, x2s - 1 - gind2, ii) = quan(specIND, valIND, kk, x2s + gind2, ii);
                    }
                }
            }
        }
    }
}

void dust_standard_boundary_condition_x2o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    #pragma omp parallel for collapse(4)
    for (int specIND = 0; specIND < quan.shape()[0]; specIND ++){
        for (int valIND = 0; valIND < quan.shape()[1]; valIND ++){
            for (int gind2 = 0; gind2 < ng2; gind2 ++){
                for (int kk = x3s; kk < x3l; kk++){
                    for (int ii = x1s; ii < x1l; ii++){
                        quan(specIND, valIND, kk, x2l + gind2, ii)     = quan(specIND, valIND, kk, x2l - (gind2 + 1), ii);
                    }
                }
            }
        }
    }
}


void dust_standard_boundary_condition_x3i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    #pragma omp parallel for collapse(4)
    for (int specIND = 0; specIND < quan.shape()[0]; specIND ++){
        for (int valIND = 0; valIND < quan.shape()[1]; valIND ++){
            for (int gind3 = 0; gind3 < ng3; gind3 ++){
                for (int jj = x2s; jj < x2l; jj++){
                    for (int ii = x1s; ii < x1l; ii++){
                        quan(specIND, valIND, x3s - 1 - gind3, jj, ii) = quan(specIND, valIND, x3s + gind3, jj, ii);
                    }
                }
            }
        }
    }
}


void dust_standard_boundary_condition_x3o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3){
    // standard inflow/outflow boundary, simply copy the values from active zone to ghost zone.
    #pragma omp parallel for collapse(4)
    for (int specIND = 0; specIND < quan.shape()[0]; specIND ++){
        for (int valIND = 0; valIND < quan.shape()[1]; valIND ++){
            for (int gind3 = 0; gind3 < ng3; gind3 ++){
                for (int jj = x2s; jj < x2l; jj++){
                    for (int ii = x1s; ii < x1l; ii++){
                        quan(specIND, specIND, valIND, x3l + gind3, jj, ii)     = quan(specIND, specIND, valIND, x3l - (gind3 + 1), jj, ii);
                    }
                }
            }
        }
    }
}

