#include "../algorithm/mesh.hpp"
#include "../algorithm/YishengArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include <cstdlib>
#include "../algorithm/inoutput/input.hpp"

void setup(mesh &m, input_file &finput){
    float gamma_hydro = 1.4;
    m.SetupCartesian(3,                           // 3D problem
                     -0.5, 0.5, 64, 2,            // ax1, index 2
                     -0.5, 0.5, 64, 2,            // ax2, index 1
                     -0.5, 0.5, 64, 2             // ax3, index 0
                     );
    m.hydro_gamma = gamma_hydro;

    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                if (pow(m.x1v(ii), 2) + pow(m.x2v(jj), 2) + pow(m.x3v(kk), 2) < pow(0.1, 2)){
                    m.cons(IDN, kk, jj, ii) = 5.;
                    m.cons(IM1, kk, jj, ii) = 0;
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    m.prim(IPN, kk, jj, ii) = 4.5;
                    float vel1 = m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel2 = m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel3 = m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    m.cons(IEN, kk, jj, ii) = ene(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii),
                                                                  vel1,
                                                                  vel2,
                                                                  vel3, m.hydro_gamma);
                }
                else{
                    m.cons(IDN, kk, jj, ii) = 2.0;
                    m.cons(IM1, kk, jj, ii) = 0;
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    m.prim(IPN, kk, jj, ii) = 2.5;
                    float vel1 = m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel2 = m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel3 = m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    m.cons(IEN, kk, jj, ii) = ene(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii),
                                                                  vel1,
                                                                  vel2,
                                                                  vel3, m.hydro_gamma);
                }
            }
        }
    }
    //for (int jj = m.x2s; jj < m.x2l; jj ++){
    //    m.mom1(m.nghost, jj, m.nghost) = 3;
    //}
    //m.mom1(1, m.x2s + m.nx2 * 1 / 4 + 2, 101) = 10;

    m.cons_to_prim();
    apply_boundary_condition(m);
}
