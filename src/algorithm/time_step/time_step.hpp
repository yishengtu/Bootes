
#ifndef TIME_STEP_HPP_
#define TIME_STEP_HPP_

#include "../eos/eos.hpp"
#include <limits>

float timestep(mesh &m, float &CFL){
    float min_dt = std::numeric_limits<float>::max();
    #pragma omp parallel for collapse(3) reduction (min : min_dt)
    for (int kk = m.x3s; kk < m.x3l ; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                float cs = soundspeed(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii), m.hydro_gamma);
                float vmx1 = abs(max(cs + m.prim(IV1, kk, jj, ii), cs - m.prim(IV1, kk, jj, ii)));
                float vmx2 = abs(max(cs + m.prim(IV2, kk, jj, ii), cs - m.prim(IV2, kk, jj, ii)));
                float vmx3 = abs(max(cs + m.prim(IV3, kk, jj, ii), cs - m.prim(IV3, kk, jj, ii)));
                min_dt = min(min(CFL * m.dx2p(jj) / vmx2, CFL * m.dx3p(kk) / vmx3), CFL * m.dx1p(ii) / vmx1);
            }
        }
    }
    return min_dt;
}

#endif // TIME_STEP_HPP_
