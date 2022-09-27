#include "const_reconst_dust.hpp"
#include "../BootesArray.hpp"
#include "../eos/momentum.hpp"
#include "../../defs.hpp"
#include "../index_def.hpp"
#include "../mesh/mesh.hpp"


void reconstruct_dust_const(mesh &m,
                            BootesArray<double> &valsL,
                            BootesArray<double> &valsR,
                            int &x1excess, int &x2excess, int &x3excess,
                            int &axis,
                            int &IMP,
                            double &dt
                            ){
    // Computation starts in first ghost zone, for first active cell left boundary flux
   // #pragma omp parallel for collapse (3) schedule (static) firstprivate(dt)
    #pragma acc parallel loop collapse(4) default(present) firstprivate(x1excess,x2excess,x3excess,axis)
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
        for (int kk = -x3excess; kk < m.nx3 + x3excess; kk++){
            for (int jj = -x2excess; jj < m.nx2 + x2excess; jj++){
                for (int ii = -x1excess; ii < m.nx1 + x1excess; ii++){
                    // Left of a cell is the right of an edge.
                    if (kk == -1 || jj == -1 || ii == -1){
                        ;
                    }
                    else {
                        valsR(specIND, axis, IDN, kk, jj, ii) = m.dcons(specIND, IDN, m.x3s + kk, m.x2s + jj, m.x1s + ii);
                        valsR(specIND, axis, IM1, kk, jj, ii) = m.dcons(specIND, IM1, m.x3s + kk, m.x2s + jj, m.x1s + ii);
                        valsR(specIND, axis, IM2, kk, jj, ii) = m.dcons(specIND, IM2, m.x3s + kk, m.x2s + jj, m.x1s + ii);
                        valsR(specIND, axis, IM3, kk, jj, ii) = m.dcons(specIND, IM3, m.x3s + kk, m.x2s + jj, m.x1s + ii);
                    }

                    if (kk == m.nx3 || jj == m.nx2 || ii == m.nx1){
                        ;
                    }
                    else{
                        valsL(specIND, axis, IDN, kk + x3excess, jj + x2excess, ii + x1excess) = m.dcons(specIND, IDN, m.x3s + kk, m.x2s + jj, m.x1s + ii);
                        valsL(specIND, axis, IM1, kk + x3excess, jj + x2excess, ii + x1excess) = m.dcons(specIND, IM1, m.x3s + kk, m.x2s + jj, m.x1s + ii);
                        valsL(specIND, axis, IM2, kk + x3excess, jj + x2excess, ii + x1excess) = m.dcons(specIND, IM2, m.x3s + kk, m.x2s + jj, m.x1s + ii);
                        valsL(specIND, axis, IM3, kk + x3excess, jj + x2excess, ii + x1excess) = m.dcons(specIND, IM3, m.x3s + kk, m.x2s + jj, m.x1s + ii);
                    }
                }
            }
        }
    }
}

