#include "gas_drag_on_dust.hpp"
#include "../mesh/mesh.hpp"
#include "../index_def.hpp"
#include "../eos/eos.hpp"


double stoppingtime(double &rhodmsize, double &rho, double &pres, double &vth_coeff){
    // cout << rhodmsize << '\t' << rhodmsize / (rho * thermalspeed(rho, pres, vth_coeff)) << endl << flush;
    return rhodmsize / (rho * thermalspeed(rho, pres, vth_coeff));
}


double calc_stoppingtimemesh(mesh &m, BootesArray<double> &stoppingtimemesh){
    #pragma omp parallel for collapse (4)
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    stoppingtimemesh(specIND, kk, jj, ii) = stoppingtime(m.GrainSizeTimesGrainDensity(specIND), m.prim(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii), m.vth_coeff);
                }
            }
        }
    }
    double minstoppingtime = std::numeric_limits<double>::max();
    #pragma omp parallel for reduction (min : minstoppingtime)
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    minstoppingtime = min(minstoppingtime, stoppingtimemesh(specIND, kk, jj, ii));
                }
            }
        }
    }
    return minstoppingtime;
}


void apply_gas_drag_on_dust(mesh &m, double &dt, BootesArray<double> &stoppingtimemesh){
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    double vdust1 = m.dprim(specIND, IV1, kk, jj, ii);
                    double vgas1  = m.prim(IV1, kk, jj, ii);
                    double vdust2 = m.dprim(specIND, IV2, kk, jj, ii);
                    double vgas2  = m.prim(IV2, kk, jj, ii);
                    double vdust3 = m.dprim(specIND, IV3, kk, jj, ii);
                    double vgas3  = m.prim(IV3, kk, jj, ii);
                    double rhodt_stime = m.dprim(specIND, IDN, kk, jj, ii) * dt / stoppingtimemesh(specIND, kk, jj, ii);
                    double dragMOM1 = rhodt_stime * (vgas1 - vdust1);
                    double dragMOM2 = rhodt_stime * (vgas2 - vdust2);
                    double dragMOM3 = rhodt_stime * (vgas3 - vdust3);
                    m.dcons(specIND, IM1, kk, jj, ii) += dragMOM1;
                    m.dcons(specIND, IM2, kk, jj, ii) += dragMOM2;
                    m.dcons(specIND, IM3, kk, jj, ii) += dragMOM3;
                }
            }
        }
    }
}


