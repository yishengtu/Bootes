#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include <cstdlib>
#include "../algorithm/inoutput/input.hpp"

void setup(mesh &m, input_file &finput){
    srand(time(0));
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                if (abs(m.x2v(jj)) < 0.1 && abs(m.x3v(kk)) < 0.1){
                    m.cons(IDN, kk, jj, ii) = 2.0;
                    m.cons(IM1, kk, jj, ii) = - 0.5 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    if (m.x1s + m.nx1 * 2 / 5 - 1 < ii && ii < m.x1s + m.nx1 * 3 / 5){
                        m.cons(IM2, kk, jj, ii) += 0.1 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                        m.cons(IM3, kk, jj, ii) += 0.1 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    }
                    m.prim(IPN, kk, jj, ii) = 2.5;
                    double vel1 = m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    double vel2 = m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    double vel3 = m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    m.cons(IEN, kk, jj, ii) = ene(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii),
                                                                  vel1,
                                                                  vel2,
                                                                  vel3, m.hydro_gamma);
                }
                else{
                    m.cons(IDN, kk, jj, ii) = 1.0;
                    m.cons(IM1, kk, jj, ii) = 0.5 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    if (m.x1s + m.nx1 * 3 / 7 - 1 < ii && ii < m.x1s + m.nx1 * 4 / 7){
                        m.cons(IM2, kk, jj, ii) += 0.1 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                        m.cons(IM3, kk, jj, ii) += 0.1 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    }
                    m.prim(IPN, kk, jj, ii) = 2.5;
                    double vel1 = m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    double vel2 = m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    double vel3 = m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    m.cons(IEN, kk, jj, ii) = ene(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii),
                                                                  vel1,
                                                                  vel2,
                                                                  vel3, m.hydro_gamma);
                }
            }
        }
    }
}


void work_after_loop(mesh &m){
    ;
}


