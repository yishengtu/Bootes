#include <limits>
#include "gas_drag_on_dust.hpp"
#include "../mesh/mesh.hpp"
#include "../index_def.hpp"
#include "../eos/eos.hpp"


double stoppingtime(double &rhodmsize, double &rho, double &pres, double &vth_coeff){
    // cout << rhodmsize << '\t' << rhodmsize / (rho * thermalspeed(rho, pres, vth_coeff)) << endl << flush;
    return rhodmsize / (rho * thermalspeed(rho, pres, vth_coeff));
}


void calc_stoppingtimemesh(mesh &m, BootesArray<double> &stoppingtimemesh){
    //#pragma omp parallel for collapse (4)
    #pragma acc parallel loop collapse (4) default (present)
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    stoppingtimemesh(specIND, kk, jj, ii) = stoppingtime(m.GrainSizeTimesGrainDensity(specIND), m.prim(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii), m.vth_coeff);
                }
            }
        }
    }
}


double find_smallest_stoppingtime(mesh &m, BootesArray<double> &stoppingtimemesh){
    double minstoppingtime = std::numeric_limits<double>::max();
    #pragma omp parallel for reduction (min : minstoppingtime)
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    minstoppingtime = std::min(minstoppingtime, stoppingtimemesh(specIND, kk, jj, ii));
                }
            }
        }
    }
    return minstoppingtime;
}
