#include <cmath>

//#include "../reconstruct/const_recon.hpp"
#include "../reconstruct/minmod.hpp"
//#include "../reconstruct/MUSCL_Hancock.hpp"
#include "../time_step/time_step.hpp"
#include "../BootesArray.hpp"
#include "../util/util.hpp"
//#include "../hydro/hll.hpp"
#include "../hydro/hlle.hpp"
//#include "../hydro/hllc.hpp"
#include "../boundary_condition/apply_bc.hpp"
#include "../index_def.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"


void calc_flux(mesh &m, double dt){
    // store the reconstructed value
    // index: (advecting direction, quantity, kk, jj, ii)
    // #pragma acc parallel loop gang default (present) firstprivate(dt)
    for (int axis = 0; axis < m.dim; axis ++){
        // step 1.1: reconstruct left/right values
        int x1excess, x2excess, x3excess;
        int IMP;        // Index of the velocity used in this axis
        if      (axis == 0){ x1excess = 1; x2excess = 0; x3excess = 0; IMP = IM1;}
        else if (axis == 1){ x1excess = 0; x2excess = 1; x3excess = 0; IMP = IM2;}
        else if (axis == 2){ x1excess = 0; x2excess = 0; x3excess = 1; IMP = IM3;}
        else { ; }

        reconstruct_minmod(m, m.valsL, m.valsR, x1excess, x2excess, x3excess, axis, IMP, dt);
        // update to device
        m.valsL.updatedev();
        m.valsR.updatedev();
        // step 1.2: solve the Riemann problem. Use HLL for now, update conservative vars
        //#pragma omp parallel for collapse (3) schedule (static)
        // present: if it's already there, don't do anything.
        // #pragma acc loop collapse (3) vector // default (present)
        #pragma acc parallel loop collapse (3) default (present)
        for (int kk = 0; kk < m.nx3 + x3excess; kk ++){
            for (int jj = 0; jj < m.nx2 + x2excess; jj ++){
                // #pragma acc loop vector private (valL, valR, fxs)   // parallel nx1 on the x-dimension of the cuda block
                for (int ii = 0; ii < m.nx1 + x1excess; ii ++){
                    double valL[5];
                    double valR[5];
                    double fxs[5];
                    valL[IDN] = m.valsL(axis, IDN, kk, jj, ii); valR[IDN] = m.valsR(axis, IDN, kk, jj, ii);
                    valL[IM1] = m.valsL(axis, IM1, kk, jj, ii); valR[IM1] = m.valsR(axis, IM1, kk, jj, ii);
                    valL[IM2] = m.valsL(axis, IM2, kk, jj, ii); valR[IM2] = m.valsR(axis, IM2, kk, jj, ii);
                    valL[IM3] = m.valsL(axis, IM3, kk, jj, ii); valR[IM3] = m.valsR(axis, IM3, kk, jj, ii);
                    valL[IEN] = m.valsL(axis, IEN, kk, jj, ii); valR[IEN] = m.valsR(axis, IEN, kk, jj, ii);
                    #ifdef ENABLE_TEMPERATURE_PROTECTION
                    valL[IEN] = energy_from_temperature_protection(valL[IDN], valL[IEN], valL[IM1], valL[IM2], valL[IM3], m.minTemp, m.hydro_gamma);
                    valR[IEN] = energy_from_temperature_protection(valR[IDN], valR[IEN], valR[IM1], valR[IM2], valR[IM3], m.minTemp, m.hydro_gamma);
                    #endif // ENABLE_TEMPERATURE_PROTECTION
                    hlle(valL, valR, fxs,
                         IMP,                   // the momentum term to add pressure; shift by one index (since first index is density)
                         m.hydro_gamma
                         );
                    m.fcons(IDN, axis, kk, jj, ii) = fxs[IDN];
                    m.fcons(IM1, axis, kk, jj, ii) = fxs[IM1];
                    m.fcons(IM2, axis, kk, jj, ii) = fxs[IM2];
                    m.fcons(IM3, axis, kk, jj, ii) = fxs[IM3];
                    m.fcons(IEN, axis, kk, jj, ii) = fxs[IEN];
                }
            }
        }
    }
    // Need to set unused values in the fdcons to zeros.
    // To do so the axis goes from "number of active axis" to 3
    // #pragma omp parallel for collapse (4) schedule (static)
    #pragma acc parallel loop collapse (4) default (present)
    for (int axis = m.dim; axis < 3; axis ++){
        for (int kk = 0; kk < m.fcons.shape()[2]; kk ++){
            for (int jj = 0; jj < m.fcons.shape()[3]; jj ++){
                for (int ii = 0; ii < m.fcons.shape()[4]; ii ++){
                    m.fcons(IDN, axis, kk, jj, ii) = 0.0;
                    m.fcons(IM1, axis, kk, jj, ii) = 0.0;
                    m.fcons(IM2, axis, kk, jj, ii) = 0.0;
                    m.fcons(IM3, axis, kk, jj, ii) = 0.0;
                    m.fcons(IEN, axis, kk, jj, ii) = 0.0;
                }
            }
        }
    }
    // copy m.fcons back to host
    m.fcons.updatehost();
}


void advect_cons(mesh &m, double dt){
    #if defined(CARTESIAN_COORD)
        //#pragma omp parallel for collapse (3) schedule (static)
        #pragma acc parallel loop collapse (3) default(present)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    int kkf = kk - m.x3s;
                    int jjf = jj - m.x2s;
                    int iif = ii - m.x1s;
                    for (int consIND = 0; consIND < NUMCONS; consIND++){
                        m.cons(consIND, kk, jj, ii) -= (dt / m.dx1(ii) * (m.fcons(consIND, 0, kkf, jjf, iif + 1) - m.fcons(consIND, 0, kkf, jjf, iif))
                                                      + dt / m.dx2(jj) * (m.fcons(consIND, 1, kkf, jjf + 1, iif) - m.fcons(consIND, 1, kkf, jjf, iif))
                                                      + dt / m.dx3(kk) * (m.fcons(consIND, 2, kkf + 1, jjf, iif) - m.fcons(consIND, 2, kkf, jjf, iif)));
                    }
                }
            }
        }
    #elif defined(SPHERICAL_POLAR_COORD)
        // #pragma omp parallel for collapse (3) schedule (static)
        #pragma acc parallel loop collapse (3) default(present)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    int kkf = kk - m.x3s;
                    int jjf = jj - m.x2s;
                    int iif = ii - m.x1s;
                    // first take care of the divergent terms
                    for (int consIND = 0; consIND < NUMCONS; consIND++){
                        m.cons(consIND, kk, jj, ii) -= (dt / m.vol(kk, jj, ii) * (m.fcons(consIND, 0, kkf, jjf, iif + 1) * m.f1a(kk, jj, ii + 1) - m.fcons(consIND, 0, kkf, jjf, iif) * m.f1a(kk, jj, ii))
                                                      + dt / m.vol(kk, jj, ii) * (m.fcons(consIND, 1, kkf, jjf + 1, iif) * m.f2a(kk, jj + 1, ii) - m.fcons(consIND, 1, kkf, jjf, iif) * m.f2a(kk, jj, ii))
                                                      + dt / m.vol(kk, jj, ii) * (m.fcons(consIND, 2, kkf + 1, jjf, iif) * m.f3a(kk + 1, jj, ii) - m.fcons(consIND, 2, kkf, jjf, iif) * m.f3a(kk, jj, ii)));
                    }
                    // geometry term
                    // lower order, but easier to understand..
                    m.cons(IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * \
                                (m.prim(IDN, kk, jj, ii) * (pow(m.prim(IV2, kk, jj, ii), 2) + pow(m.prim(IV3, kk, jj, ii), 2)) + 2 * m.prim(IPN, kk, jj, ii));
                    m.cons(IM2, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.prim(IDN, kk, jj, ii) * m.prim(IV1, kk, jj, ii) * m.prim(IV2, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) += dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * pow(m.prim(IV3, kk, jj, ii), 2) + m.prim(IPN, kk, jj, ii));
                    m.cons(IM3, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.prim(IDN, kk, jj, ii) * m.prim(IV1, kk, jj, ii) * m.prim(IV3, kk, jj, ii);
                    m.cons(IM3, kk, jj, ii) -= dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * (m.prim(IV2, kk, jj, ii) * m.prim(IV3, kk, jj, ii)));

                    /*
                    m.cons(IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * (pow(m.prim(IV2, kk, jj, ii), 2) + pow(m.prim(IV3, kk, jj, ii), 2)) + 2 * m.prim(IPN, kk, jj, ii));

                    m.cons(IM2, kk, jj, ii) -= dt * m.dx1(ii) / m.rV(ii) * (m.rsq(ii) * m.valsL(0, IM2, kkf, jjf, iif) + m.rsq(ii + 1) * m.valsR(0, IM2, kkf, jjf, iif + 1));
                    m.cons(IM2, kk, jj, ii) += dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * pow(m.prim(IV3, kk, jj, ii), 2) + m.prim(IPN, kk, jj, ii));

                    m.cons(IM3, kk, jj, ii) -= dt * m.dx1(ii) / m.rV(ii) * (m.rsq(ii) * m.valsL(0, IM3, kkf, jjf, iif) + m.rsq(ii + 1) * m.valsR(0, IM3, kkf, jjf, iif + 1));
                    m.cons(IM3, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.geo_cot(jj) / (m.geo_sm(jj) + m.geo_sp(jj)) *
                                                            (m.geo_sm(jj) * m.valsL(1, IM2, kkf, jjf, iif) * m.valsL(1, IM3, kkf, jjf, iif)
                                                            + m.geo_sp(jj) * m.valsR(1, IM2, kkf, jjf + 1, iif) * m.valsR(1, IM3, kkf, jjf + 1, iif));
                    */
                }
            }
        }
    #else
        # error need coordinate defined
    #endif
}


#ifdef DENSITY_PROTECTION
void protection(mesh &m){
    //#pragma omp parallel for collapse (3)
    #pragma acc parallel loop collapse (3) default(present)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                if (isnan(m.cons(IDN, kk, jj, ii)) || m.cons(IDN, kk, jj, ii) < m.minDensity){
                    m.cons(IDN, kk, jj, ii) = m.minDensity;
                }
                if (isnan(m.cons(IM1, kk, jj, ii))){
                    m.cons(IM1, kk, jj, ii) = m.minDensity * m.prim(IM1, kk, jj, ii);
                }
                if (isnan(m.cons(IM2, kk, jj, ii))){
                    m.cons(IM2, kk, jj, ii) = m.minDensity * m.prim(IM2, kk, jj, ii);
                }
                if (isnan(m.cons(IM3, kk, jj, ii))){
                    m.cons(IM3, kk, jj, ii) = m.minDensity * m.prim(IM3, kk, jj, ii);
                }
                if (isnan(m.cons(IEN, kk, jj, ii)) || m.cons(IEN, kk, jj, ii) < m.minDensity){
                    m.cons(IEN, kk, jj, ii) = m.minDensity;
                }
            }
        }
    }
}
#endif // DENSITY_PROTECTION

#ifdef ENABLE_TEMPERATURE_PROTECTION
void temperature_protection(mesh &m){
    //#pragma omp parallel for collapse (3)
    #pragma acc parallel loop collapse (3) default(present)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                m.cons(IEN, kk, jj, ii) = energy_from_temperature_protection(m.cons(IDN, kk, jj, ii), m.cons(IEN, kk, jj, ii),
                                                                             m.cons(IM1, kk, jj, ii), m.cons(IM2, kk, jj, ii), m.cons(IM3, kk, jj, ii), m.minTemp, m.hydro_gamma);
            }
        }
    }
}
#endif // ENABLE_TEMPERATURE_PROTECTION


