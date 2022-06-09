#include <cstdlib>

#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include "../algorithm/inoutput/input.hpp"
#include "../algorithm/gravity/gravity.hpp"
#include "../algorithm/index_def.hpp"

void setup(mesh &m, input_file &finput){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                // if (abs(m.x1v(ii) * cos(m.x2v(jj))) < 0.3){
                if (abs(((double) jj) - ((double) m.x2v.shape()[0] - 1) / 2.) < ((double) m.x2v.shape()[0]) / 16. && abs(m.x1v(ii) - 6.0) < 0.1){
                    m.cons(IDN, kk, jj, ii) = 1.0;
                    m.cons(IM1, kk, jj, ii) = 0.0; // 5. * m.cons(IDN, kk, jj, ii) * sin(m.x2v(jj));
                    m.cons(IM2, kk, jj, ii) = 0.0; // 5. * m.cons(IDN, kk, jj, ii) * cos(m.x2v(jj));
                    //m.cons(IM1, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    //m.cons(IM2, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM3, kk, jj, ii) = 0.0; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                    m.prim(IPN, kk, jj, ii) = 5.5;
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
                    m.cons(IM3, kk, jj, ii) = 0.0; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
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

    for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    m.dcons(specIND, IDN, kk, jj, ii) = 1.; // (double) specIND + 1.;
                    m.dcons(specIND, IM1, kk, jj, ii) = 0.0;
                    m.dcons(specIND, IM2, kk, jj, ii) = 0.0;
                    m.dcons(specIND, IM3, kk, jj, ii) = 0.0; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                    m.dprim(specIND, IDN, kk, jj, ii) = 1.; // (double) specIND + 1.;
                    m.dprim(specIND, IM1, kk, jj, ii) = 0.0;
                    m.dprim(specIND, IM2, kk, jj, ii) = 0.0;
                    m.dprim(specIND, IM3, kk, jj, ii) = 0.0; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                }
            }
        }
    }

    /** gravity **/
    // Right now, gravity is defined in main.cpp and time_integration.cpp.
}


void work_after_loop(mesh &m, double &dt){
    /*
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                // if (abs(m.x1v(ii) * cos(m.x2v(jj))) < 0.3){
                if (abs(((double) jj) - ((double) m.x2v.shape()[0] - 1) / 2.) < ((double) m.x2v.shape()[0]) / 16. && abs(m.x1v(ii) - 1.0) < 0.1){
                    m.cons(IDN, kk, jj, ii) = 2.0;
                    m.cons(IM1, kk, jj, ii) = 5. * m.cons(IDN, kk, jj, ii) * sin(m.x2v(jj));
                    m.cons(IM2, kk, jj, ii) = 5. * m.cons(IDN, kk, jj, ii) * cos(m.x2v(jj));
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
                    for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
                        m.dcons(specIND, IDN, kk, jj, ii) = 1.; // (double) specIND + 1.;
                    }
                }
            }
        }
    }
    */
}

void setup_dust(mesh &m, input_file &finput){
    double smin    = finput.getDouble("smin");
    double smax    = finput.getDouble("smax");
    double rhodm   = finput.getDouble("srho");
    int ns         = finput.getInt("num_species");
    m.NUMSPECIES = ns;
    m.rhodm = rhodm;

    m.GrainEdgeList = logspace(log10(smin), log10(smax), ns + 1, true);
    m.GrainSizeList.NewBootesArray(ns);
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
        double s2 = m.GrainEdgeList(specIND + 1);
        double s1 = m.GrainEdgeList(specIND);
        m.GrainSizeList(specIND) = pow((pow(s2, 4) - pow(s1, 4)) / (4 * (s2 - s1)), 1./3.);
    }
}
