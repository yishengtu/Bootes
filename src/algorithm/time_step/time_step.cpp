#include "time_step.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../index_def.hpp"
#include <limits>

double timestep(mesh &m, double &CFL){
    double min_dt = std::numeric_limits<double>::max();
    if (m.dim == 1){
        #pragma omp parallel for collapse(3) reduction (min : min_dt)
        for (int kk = m.x3s; kk < m.x3l ; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double cs = soundspeed(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii), m.hydro_gamma);
                    double vmx1 = abs(std::max(cs + m.prim(IV1, kk, jj, ii), cs - m.prim(IV1, kk, jj, ii)));

                    double dx1_sig = std::min(std::min(m.dx1p(kk, jj, ii - 1), m.dx1p(kk, jj, ii)), m.dx1p(kk, jj, ii + 1));
                    min_dt = std::min(dx1_sig / vmx1, min_dt);
                }
            }
        }
    }
    else if (m.dim == 2){
        #pragma omp parallel for collapse(3) reduction (min : min_dt)
        for (int kk = m.x3s; kk < m.x3l ; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double cs = soundspeed(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii), m.hydro_gamma);
                    double vmx1 = abs(std::max(cs + m.prim(IV1, kk, jj, ii), cs - m.prim(IV1, kk, jj, ii)));
                    double vmx2 = abs(std::max(cs + m.prim(IV2, kk, jj, ii), cs - m.prim(IV2, kk, jj, ii)));

                    double dx1_sig = std::min(std::min(m.dx1p(kk, jj, ii - 1), m.dx1p(kk, jj, ii)), m.dx1p(kk, jj, ii + 1));
                    double dx2_sig = std::min(std::min(m.dx2p(kk, jj - 1, ii), m.dx2p(kk, jj, ii)), m.dx2p(kk, jj + 1, ii));
                    min_dt = std::min(std::min(dx2_sig / vmx2, dx1_sig / vmx1), min_dt);
                    if (std::min(dx2_sig / vmx2, dx1_sig / vmx1) < 0 || min_dt < 1e-100){
                        std::cout << m.cons(IDN, kk, jj, ii) << '\t' << m.prim(IPN, kk, jj, ii) << '\t' <<  m.hydro_gamma << '\t' << cs << '\t' << m.prim(IV1, kk, jj, ii) << std::endl << std::flush;
                        std::cout << kk << '\t' << jj << '\t' << ii << '\t' << dx2_sig << '\t' << vmx2 << '\t' << dx1_sig << '\t' << vmx1 << std::endl << std::flush;
                        throw 1;
                    }
                }
            }
        }
        /*
        #ifdef ENABLE_DUSTFLUID
        #pragma omp parallel for collapse(4) reduction (min: min_dt)
        for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
            for (int kk = m.x3s; kk < m.x3l ; kk++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    for (int ii = m.x1s; ii < m.x1l; ii++){
                        double vmx1 = abs(m.dprim(specIND, IV1, kk, jj, ii));
                        double vmx2 = abs(m.dprim(specIND, IV2, kk, jj, ii));

                        double dx1_sig = std::min(std::min(m.dx1p(kk, jj, ii - 1), m.dx1p(kk, jj, ii)), m.dx1p(kk, jj, ii + 1));
                        double dx2_sig = std::min(std::min(m.dx2p(kk, jj - 1, ii), m.dx2p(kk, jj, ii)), m.dx2p(kk, jj + 1, ii));
                        min_dt = std::min(std::min(dx2_sig / vmx2, dx1_sig / vmx1), min_dt);
                    }
                }
            }
        }
        #endif // ENABLE_DUSTFLUID
        */
    }
    else{
        #pragma omp parallel for collapse(3) reduction (min : min_dt)
        for (int kk = m.x3s; kk < m.x3l ; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double cs = soundspeed(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii), m.hydro_gamma);
                    double vmx1 = abs(std::max(cs + m.prim(IV1, kk, jj, ii), cs - m.prim(IV1, kk, jj, ii)));
                    double vmx2 = abs(std::max(cs + m.prim(IV2, kk, jj, ii), cs - m.prim(IV2, kk, jj, ii)));
                    double vmx3 = abs(std::max(cs + m.prim(IV3, kk, jj, ii), cs - m.prim(IV3, kk, jj, ii)));

                    double dx1_sig = std::min(std::min(m.dx1p(kk, jj, ii - 1), m.dx1p(kk, jj, ii)), m.dx1p(kk, jj, ii + 1));
                    double dx2_sig = std::min(std::min(m.dx2p(kk, jj - 1, ii), m.dx2p(kk, jj, ii)), m.dx2p(kk, jj + 1, ii));
                    double dx3_sig = std::min(std::min(m.dx3p(kk - 1, jj, ii), m.dx3p(kk, jj, ii)), m.dx3p(kk + 1, jj, ii));
                    min_dt = std::min(std::min(std::min(dx2_sig / vmx2, dx3_sig / vmx3), dx1_sig / vmx1), min_dt);
                }
            }
        }
    }

    return CFL * min_dt;
}

