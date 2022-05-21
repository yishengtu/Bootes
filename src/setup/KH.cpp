#include "../algorithm/mesh.hpp"
#include "../algorithm/YishengArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include <cstdlib>
#include "../algorithm/inoutput/input.hpp"


void setup(mesh &m, input_file &finput){
    float gamma_hydro = 1.4;
    m.SetupCartesian(2,
                     -2.,  2., 256, 2,                       // ax1
                     -1.,  1., 128, 2,                       // ax2
                      -1.,  1.,  1,  0                        // ax3
                     );
    m.hydro_gamma = gamma_hydro;
    srand(time(0));
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                if (abs(m.x2v(jj)) < 0.1){
                    m.cons(IDN, kk, jj, ii) = 2.0;
                    m.cons(IM1, kk, jj, ii) = - 0.5 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM1, kk, jj, ii) += 0.01 * ((float) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) += 0.01 * ((float) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
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
                else{
                    m.cons(IDN, kk, jj, ii) = 1.0;
                    m.cons(IM1, kk, jj, ii) = 0.5 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM1, kk, jj, ii) += 0.01 * ((float) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) += 0.01 * ((float) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
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

    m.cons_to_prim();
    apply_boundary_condition(m);
}
