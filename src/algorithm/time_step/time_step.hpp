
#ifndef TIME_STEP_HPP_
#define TIME_STEP_HPP_

#include "../eos/eos.hpp"
#include <limits>

double timestep(mesh &m, double &CFL){
    double min_dt = std::numeric_limits<double>::max();
    #pragma omp parallel for collapse(3) reduction (min : min_dt)
    for (int kk = m.x3s; kk < m.x3l ; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                double cs = soundspeed(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii), m.hydro_gamma);
                double vmx1 = abs(max(cs + m.prim(IV1, kk, jj, ii), cs - m.prim(IV1, kk, jj, ii)));
                double vmx2 = abs(max(cs + m.prim(IV2, kk, jj, ii), cs - m.prim(IV2, kk, jj, ii)));
                double vmx3 = abs(max(cs + m.prim(IV3, kk, jj, ii), cs - m.prim(IV3, kk, jj, ii)));

                double dx1_sig = min(min(m.dx1p(ii - 1), m.dx1p(ii)), m.dx1p(ii + 1));
                double dx2_sig = min(min(m.dx2p(ii - 1), m.dx2p(ii)), m.dx2p(ii + 1));
                double dx3_sig = min(min(m.dx3p(ii - 1), m.dx3p(ii)), m.dx3p(ii + 1));
                min_dt = min(min(CFL * dx2_sig / vmx2, CFL * dx3_sig / vmx3), CFL * dx1_sig / vmx1);
            }
        }
    }
    return min_dt;
}

#endif // TIME_STEP_HPP_
