#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include <cstdlib>
#include "../algorithm/inoutput/input.hpp"


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


void setup_dust_IC(mesh &m, input_file &finput){
    double amin = finput.getDouble("ainimin");      // m.GrainEdgeList(0);
    double amax = finput.getDouble("ainimax");      // m.GrainEdgeList(m.NUMSPECIES);
    m.dminDensity = finput.getDouble("dminDensity");

    double Adust = 3. / (100. * 8 * M_PI * m.rhodm * (sqrt(amax) - sqrt(amin)));

    for (int ss = 0; ss < m.NUMSPECIES; ss++){
        if (m.GrainEdgeList(ss + 1) < amin || m.GrainEdgeList(ss) > amax){
            for (int kk = m.x3s; kk < m.x3l; kk++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    for (int ii = m.x1s; ii < m.x1l; ii++){
                        m.dcons(ss, IDN, kk, jj, ii) = m.dminDensity;
                        m.dcons(ss, IM1, kk, jj, ii) = m.dminDensity * m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        m.dcons(ss, IM2, kk, jj, ii) = m.dminDensity * m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        m.dcons(ss, IM3, kk, jj, ii) = m.dminDensity * m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    }
                }
            }
        }
        else{
            for (int kk = m.x3s; kk < m.x3l; kk++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    for (int ii = m.x1s; ii < m.x1l; ii++){
                        // double upperlim = min(amax, m.GrainEdgeList(ss + 1));
                        // double lowerlim = max(amin, m.GrainEdgeList(ss));
                        double rho = 1.0; // Adust * m.cons(IDN, kk, jj, ii) * (sqrt(upperlim) - sqrt(lowerlim)) / (0.5 * pow(m.GrainSizeList(ss), 3)) * m.GrainMassList(ss);
                        m.dcons(ss, IDN, kk, jj, ii) = 1.0;
                        m.dcons(ss, IM1, kk, jj, ii) = rho * m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        m.dcons(ss, IM2, kk, jj, ii) = rho * m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        m.dcons(ss, IM3, kk, jj, ii) = rho * m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    }
                }
            }
        }
    }
}


void setup(mesh &m, input_file &finput){
    srand(time(0));
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                if (abs(m.x2v(jj)) < 0.5){
                    m.cons(IDN, kk, jj, ii) = 2.0;
                    m.cons(IM1, kk, jj, ii) = - 0.5 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
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
                else{
                    m.cons(IDN, kk, jj, ii) = 1.0;
                    m.cons(IM1, kk, jj, ii) = 0.5 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
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

    setup_dust_IC(m, finput);
}


void work_after_loop(mesh &m, double &dt){
    ;
}
