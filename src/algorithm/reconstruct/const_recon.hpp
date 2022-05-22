#ifndef RECONSTRUCT_HPP_
#define RECONSTRUCT_HPP_

#include "../YishengArray.hpp"

void reconstruct(mesh &m,
                   YishengArray<double> &rhoL, YishengArray<double> &rhoR,
                   YishengArray<double> &v1L, YishengArray<double> &v1R,
                   YishengArray<double> &v2L, YishengArray<double> &v2R,
                   YishengArray<double> &v3L, YishengArray<double> &v3R,
                   YishengArray<double> &pressL, YishengArray<double> &pressR,
                   YishengArray<double> &m1L, YishengArray<double> &m1R,
                   YishengArray<double> &m2L, YishengArray<double> &m2R,
                   YishengArray<double> &m3L, YishengArray<double> &m3R,
                   YishengArray<double> &eneL, YishengArray<double> &eneR,
                   int &axis,
                   double dt     // not used, just here to simplify shifting between constructions
                   ){
    int x1excess, x2excess, x3excess;
    if      (axis == 0){ x1excess = 1; x2excess = 0; x3excess = 0; }
    else if (axis == 1){ x1excess = 0; x2excess = 1; x3excess = 0; }
    else if (axis == 2){ x1excess = 0; x2excess = 0; x3excess = 1; }
    else { cout << "axis > 3!!!" << endl << flush; throw 1; }

    rhoL.NewYishengArray(  m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    rhoR.NewYishengArray(  m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    v1L.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    v1R.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    v2L.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    v2R.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    v3L.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    v3R.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    pressL.NewYishengArray(m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    pressR.NewYishengArray(m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    m1L.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    m1R.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    m2L.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    m2R.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    m3L.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    m3R.NewYishengArray(   m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    eneL.NewYishengArray(  m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);
    eneR.NewYishengArray(  m.nx3 + x3excess, m.nx2 + x2excess, m.nx1 + x1excess);

    int x1shiftL = 0;
    int x2shiftL = 0;
    int x3shiftL = 0;
    if (axis == 0){      x2shiftL = m.ng2 - 1; x3shiftL = m.ng3 - 1; }
    else if (axis == 1){ x1shiftL = m.ng1 - 1; x3shiftL = m.ng3 - 1; }
    else if (axis == 2){ x1shiftL = m.ng1 - 1; x2shiftL = m.ng2 - 1; }
    int x1shiftR = x1shiftL + x1excess;
    int x2shiftR = x2shiftL + x2excess;
    int x3shiftR = x3shiftL + x3excess;
    #pragma omp parallel for collapse (3) schedule (static)
    for (int kk = 0; kk < rhoL.shape()[0] ; kk++){
        for (int jj = 0; jj < rhoL.shape()[1]; jj++){
            for (int ii = 0; ii < rhoL.shape()[2]; ii++){
                rhoL(kk, jj, ii)   = m.rho(     kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);
                v1L(kk, jj, ii)    = m.v1(      kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);
                v2L(kk, jj, ii)    = m.v2(      kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);
                v3L(kk, jj, ii)    = m.v3(      kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);
                pressL(kk, jj, ii) = m.pressure(kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);
                m1L(kk, jj, ii)    = m.mom1(    kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);
                m2L(kk, jj, ii)    = m.mom2(    kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);
                m3L(kk, jj, ii)    = m.mom3(    kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);
                eneL(kk, jj, ii)   = m.energy(  kk + x3shiftL, jj + x2shiftL, ii + x1shiftL);

                rhoR(kk, jj, ii)   = m.rho(     kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
                v1R(kk, jj, ii)    = m.v1(      kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
                v2R(kk, jj, ii)    = m.v2(      kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
                v3R(kk, jj, ii)    = m.v3(      kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
                pressR(kk, jj, ii) = m.pressure(kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
                m1R(kk, jj, ii)    = m.mom1(    kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
                m2R(kk, jj, ii)    = m.mom2(    kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
                m3R(kk, jj, ii)    = m.mom3(    kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
                eneR(kk, jj, ii)   = m.energy(  kk + x3shiftR, jj + x2shiftR, ii + x1shiftR);
            }
        }
    }
}

#endif
