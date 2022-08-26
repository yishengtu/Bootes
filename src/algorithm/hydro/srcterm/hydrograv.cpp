#include "hydrograv.hpp"

#include "../../reconstruct/minmod.hpp"
//#include "reconstruct/const_recon.hpp"
#include "../../BootesArray.hpp"
#include "../../util/util.hpp"
#include "../hll.hpp"
#include "../../index_def.hpp"
#include "../../mesh/mesh.hpp"
#include "../../gravity/gravity.hpp"


#ifdef ENABLE_GRAVITY
void apply_grav_source_terms(mesh &m, double &dt){
    #if defined (CARTESIAN_COORD)
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    double rhogradphix1 = m.prim(IDN, kk, jj, ii) * m.grav->grav_x1(kk, jj, ii); // (m.grav->Phi_grav_x1surface(kk, jj, ii + 1) - m.grav->Phi_grav_x1surface(kk, jj, ii)) / m.dx1p(kk, jj, ii);
                    double rhogradphix2 = m.prim(IDN, kk, jj, ii) * m.grav->grav_x2(kk, jj, ii); // (m.grav->Phi_grav_x2surface(kk, jj + 1, ii) - m.grav->Phi_grav_x2surface(kk, jj, ii)) / m.dx2p(kk, jj, ii);
                    double rhogradphix3 = m.prim(IDN, kk, jj, ii) * m.grav->grav_x3(kk, jj, ii); // (m.grav->Phi_grav_x3surface(kk + 1, jj, ii) - m.grav->Phi_grav_x3surface(kk, jj, ii)) / m.dx3p(kk, jj, ii);
                    m.cons(IM1, kk, jj, ii) += rhogradphix1 * dt;
                    m.cons(IM2, kk, jj, ii) += rhogradphix2 * dt;
                    m.cons(IM3, kk, jj, ii) += rhogradphix3 * dt;
                    m.cons(IEN, kk, jj, ii) -= (rhogradphix1 * m.prim(IV1, kk, jj, ii) + rhogradphix2 * m.prim(IV2, kk, jj, ii) + rhogradphix3 * m.prim(IV3, kk, jj, ii)) * dt;
                }
            }
        }
        // TODO
    #elif defined (SPHERICAL_POLAR_COORD)
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    double rhogradphix1 = m.prim(IDN, kk, jj, ii) * m.grav->grav_x1(kk, jj, ii);
                    double rhogradphix2 = m.prim(IDN, kk, jj, ii) * m.grav->grav_x2(kk, jj, ii);
                    double rhogradphix3 = m.prim(IDN, kk, jj, ii) * m.grav->grav_x3(kk, jj, ii);
                    m.cons(IM1, kk, jj, ii) += rhogradphix1 * dt;
                    m.cons(IM2, kk, jj, ii) += rhogradphix2 * dt;
                    m.cons(IM3, kk, jj, ii) += rhogradphix3 * dt;
                    m.cons(IEN, kk, jj, ii) -= (rhogradphix1 * m.prim(IV1, kk, jj, ii) + rhogradphix2 * m.prim(IV2, kk, jj, ii) + rhogradphix3 * m.prim(IV3, kk, jj, ii)) * dt;
                }
            }
        }
    #endif // defined (coordinate)
}
#endif // ENABLE_GRAVITY
