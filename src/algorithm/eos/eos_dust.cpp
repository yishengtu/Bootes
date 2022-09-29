#include <cmath>
#include "eos_dust.hpp"
#include "../../defs.hpp"
#include "../index_def.hpp"
#include "../mesh/mesh.hpp"


#ifdef ENABLE_DUSTFLUID
void cons_to_prim_dust(mesh &m){
    #ifdef GPU
    #pragma acc parallel loop collapse (4)
    #else
    #pragma omp parallel for collapse (4) schedule (static)
    #endif
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
        for (int kk = m.x3s; kk < m.x3l ; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    m.dprim(specIND, IDN, kk, jj, ii) = m.dcons(specIND, IDN, kk, jj, ii);
                    m.dprim(specIND, IV1, kk, jj, ii) = m.dcons(specIND, IM1, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                    m.dprim(specIND, IV2, kk, jj, ii) = m.dcons(specIND, IM2, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                    m.dprim(specIND, IV3, kk, jj, ii) = m.dcons(specIND, IM3, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                }
            }
        }
    }
}

#endif // ENABLE_DUSTFLUID
