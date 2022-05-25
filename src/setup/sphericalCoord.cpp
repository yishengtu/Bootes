#include <cstdlib>

#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include "../algorithm/inoutput/input.hpp"
#include "../algorithm/gravity/gravity.hpp"

void setup(mesh &m, input_file &finput){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                // if (abs(m.x1v(ii) * cos(m.x2v(jj))) < 0.3){
                if (abs(m.x2v(jj) - M_PI / 2.) < M_PI / 16. && abs(m.x1v(ii) - 1.0) < 0.1){
                    m.cons(IDN, kk, jj, ii) = 1.0; // 2.0;
                    m.cons(IM1, kk, jj, ii) = 0.0; // 5. * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0.0; //- 0.5 * cos(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                    //m.cons(IM1, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    //m.cons(IM2, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM3, kk, jj, ii) = sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
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
                    m.cons(IM1, kk, jj, ii) = 0.0; // 0.5 * sin(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0.0; // 0.5 * cos(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                    //m.cons(IM1, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    //m.cons(IM2, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM3, kk, jj, ii) = sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
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
    gravity grav(m);
    grav.pointsource_grav(10, 0, 0, 0);
}

void work_after_loop(mesh &m){
    ;
    /*
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                // if (abs(m.x1v(ii) * cos(m.x2v(jj))) < 0.3){
                if (abs(m.x2v(jj) - M_PI / 2.) < M_PI / 20. && abs(m.x1v(ii) - 1.0) < 0.1){
                    m.cons(IDN, kk, jj, ii) = 2.0;
                    m.cons(IM1, kk, jj, ii) = 5. * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0.0; //- 0.5 * cos(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM1, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM3, kk, jj, ii) = 0;
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
    */
}
