#ifndef MUSCL_RECONSTRUCT_HPP_
#define MUSCL_RECONSTRUCT_HPP_

#include "../YishengArray.hpp"


void MH(YishengArray<double> &quan, int &x3s, int &x2s, int &x1s, double &dx_axis,
          int &kk, int &jj, int &ii, int &x3excess, int &x2excess, int &x1excess,
          int &axis, double &BquanL, double &BquanR){
    double w = 0.0;
    double Dim = quan(x3s + kk, x2s + jj, x1s + ii) - quan(x3s + kk - x3excess, x2s + jj - x2excess, x1s + ii - x1excess); // Toro 13.28
    double Dip = quan(x3s + kk + x3excess, x2s + jj + x2excess, x1s + ii + x1excess) - quan(x3s + kk, x2s + jj, x1s + ii);
    double Di = 0.5 * (1. + w) * Dim + 0.5 * (1. - w) * Dip;                   // Toro 13.27
    double rhoiL = quan(x3s + kk, x2s + jj, x1s + ii) - 0.5 * Di;         // Toro 14.33
    double rhoiR = quan(x3s + kk, x2s + jj, x1s + ii) + 0.5 * Di;         // Toro 14.33
    // In cartesian coordinate FiL = FiR because the area is the same,
    // so the flux (Toro 14.2) is the same for left and right cell boundaries
    BquanL = rhoiL; //+ 0.5 * dt / dx_axis * (FUiL - FUiR);                // Toro 14.43
    BquanR = rhoiR; //+ 0.5 * dt / dx_axis * (FUiL - FUiR);                // Toro 14.43
}


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
                   double &dt
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

    // Computation starts in first ghost zone, for first active cell left boundary flux
    #pragma omp parallel for collapse (3) schedule (static)
    for (int kk = -x3excess; kk < m.nx3 + x3excess; kk++){
        for (int jj = -x2excess; jj < m.nx2 + x2excess; jj++){
            for (int ii = -x1excess; ii < m.nx1 + x1excess; ii++){
                double dx_axis;
                if      (axis == 0) { dx_axis = m.dx3(kk); }
                else if (axis == 1) { dx_axis = m.dx2(jj); }
                else                { dx_axis = m.dx1(ii); }
                double BrhoL, BrhoR;     MH(m.rho, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, BrhoL, BrhoR);
                double Bv1L, Bv1R;       MH(m.v1, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, Bv1L, Bv1R);
                double Bv2L, Bv2R;       MH(m.v2, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, Bv2L, Bv2R);
                double Bv3L, Bv3R;       MH(m.v3, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, Bv3L, Bv3R);
                double BpressL, BpressR; MH(m.pressure, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, BpressL, BpressR);
                double Bm1L, Bm1R;       MH(m.mom1, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, Bm1L, Bm1R);
                double Bm2L, Bm2R;       MH(m.mom2, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, Bm2L, Bm2R);
                double Bm3L, Bm3R;       MH(m.mom3, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, Bm3L, Bm3R);
                double BeneL, BeneR;     MH(m.energy, m.x3s, m.x2s, m.x1s, dx_axis, kk, jj, ii, x3excess, x2excess, x1excess, axis, BeneL, BeneR);
                if (kk == -1 || jj == -1 || ii == -1){
                    ;
                }
                else {
                    rhoR(kk, jj, ii)   = BrhoL;
                    v1R(kk, jj, ii)    = Bv1L;
                    v2R(kk, jj, ii)    = Bv2L;
                    v3R(kk, jj, ii)    = Bv3L;
                    pressR(kk, jj, ii) = BpressL;
                    m1R(kk, jj, ii)    = Bm1L;
                    m2R(kk, jj, ii)    = Bm2L;
                    m3R(kk, jj, ii)    = Bm3L;
                    eneR(kk, jj, ii)   = BeneL;
                }

                if (kk == m.nx3 || jj == m.nx2 || ii == m.nx1){
                    ;
                }
                else{
                    rhoL(kk + x3excess, jj + x2excess, ii + x1excess)   = BrhoR;
                    v1L(kk + x3excess, jj + x2excess, ii + x1excess)    = Bv1R;
                    v2L(kk + x3excess, jj + x2excess, ii + x1excess)    = Bv2R;
                    v3L(kk + x3excess, jj + x2excess, ii + x1excess)    = Bv3R;
                    pressL(kk + x3excess, jj + x2excess, ii + x1excess) = BpressR;
                    m1L(kk + x3excess, jj + x2excess, ii + x1excess)    = Bm1R;
                    m2L(kk + x3excess, jj + x2excess, ii + x1excess)    = Bm2R;
                    m3L(kk + x3excess, jj + x2excess, ii + x1excess)    = Bm3R;
                    eneL(kk + x3excess, jj + x2excess, ii + x1excess)   = BeneR;
                }
            }
        }
    }
}


#endif
