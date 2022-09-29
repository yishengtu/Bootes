#include <cmath>

#include "adv_dust.hpp"

// #include "../reconstruct/minmod_dust.hpp"
#include "../reconstruct/const_reconst_dust.hpp"
#include "../time_step/time_step.hpp"
#include "../util/util.hpp"
// #include "../dust/hll_dust.hpp"
#include "../dust/doner_dust.hpp"
#include "../boundary_condition/apply_bc.hpp"
#include "../mesh/mesh.hpp"
#include "../dust/terminalvel.hpp"

void calc_flux_dust(mesh &m, double dt){
    // store the redconstructed value
    // index: (specIadvecting direction, quantity, kk, jj, ii)
    for (int axis = 0; axis < m.dim; axis ++){
        // step 1.1: redconstruct left/right values
        int x1excess, x2excess, x3excess;
        int IMP;        // Index of the velocity used in this axis
        if      (axis == 0){ x1excess = 1; x2excess = 0; x3excess = 0; IMP = IM1;}
        else if (axis == 1){ x1excess = 0; x2excess = 1; x3excess = 0; IMP = IM2;}
        else if (axis == 2){ x1excess = 0; x2excess = 0; x3excess = 1; IMP = IM3;}
        else { cout << "axis > 3!!!" << endl << flush; throw 1; }

        reconstruct_dust_const(m, m.dvalsL, m.dvalsR, x1excess, x2excess, x3excess, axis, IMP, dt);
        // step 1.2: solve the Riemann problem. Use HLL for now, update dconservative vars
        #ifdef GPU
        #pragma acc parallel loop collapse (4) default(present)
        #else
        #pragma omp parallel for collapse (4) schedule (static)
        #endif
        for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
            for (int kk = 0; kk < m.nx3 + x3excess; kk ++){
                for (int jj = 0; jj < m.nx2 + x2excess; jj ++){
                    for (int ii = 0; ii < m.nx1 + x1excess; ii ++){
                        double valL[4];
                        double valR[4];
                        double fxs[4];
                        valL[IDN] = m.dvalsL(specIND, axis, IDN, kk, jj, ii); valR[IDN] = m.dvalsR(specIND, axis, IDN, kk, jj, ii);
                        valL[IM1] = m.dvalsL(specIND, axis, IM1, kk, jj, ii); valR[IM1] = m.dvalsR(specIND, axis, IM1, kk, jj, ii);
                        valL[IM2] = m.dvalsL(specIND, axis, IM2, kk, jj, ii); valR[IM2] = m.dvalsR(specIND, axis, IM2, kk, jj, ii);
                        valL[IM3] = m.dvalsL(specIND, axis, IM3, kk, jj, ii); valR[IM3] = m.dvalsR(specIND, axis, IM3, kk, jj, ii);
                        doner_cell_dust( valL, valR, fxs,
                                         IMP,                   // the momentum term to add pressure; shift by one index (since first index is density)
                                         m.hydro_gamma
                                         );
                        m.fdcons(specIND, IDN, axis, kk, jj, ii) = fxs[IDN];
                        m.fdcons(specIND, IM1, axis, kk, jj, ii) = fxs[IM1];
                        m.fdcons(specIND, IM2, axis, kk, jj, ii) = fxs[IM2];
                        m.fdcons(specIND, IM3, axis, kk, jj, ii) = fxs[IM3];
                    }
                }
            }
        }
    }
    // Need to set unused values in the m.fdcons to zeros.
    // To do so the axis goes from "number of active axis" to 3
    #ifdef GPU
    #pragma acc parallel loop collapse (5) default(present)
    #else
    #pragma omp parallel for collapse (5) schedule (static)
    #endif
    for (int axis = m.dim; axis < 3; axis ++){
        for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
            for (int kk = 0; kk < m.fdcons.shape()[3]; kk ++){
                for (int jj = 0; jj < m.fdcons.shape()[4]; jj ++){
                    for (int ii = 0; ii < m.fdcons.shape()[5]; ii ++){
                        m.fdcons(specIND, IDN, axis, kk, jj, ii) = 0.0;
                        m.fdcons(specIND, IM1, axis, kk, jj, ii) = 0.0;
                        m.fdcons(specIND, IM2, axis, kk, jj, ii) = 0.0;
                        m.fdcons(specIND, IM3, axis, kk, jj, ii) = 0.0;
                    }
                }
            }
        }
    }
}

#if defined(SPHERICAL_POLAR_COORD)
void advect_cons_dust_sphericalpolar(mesh &m, double dt){
    #ifdef GPU
    #pragma acc parallel loop collapse (4) default(present)
    #else
    #pragma omp parallel for collapse (4) schedule (static)
    #endif
    for (int specIND  = 0; specIND < m.NUMSPECIES; specIND++){
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    int kkf = kk - m.x3s;
                    int jjf = jj - m.x2s;
                    int iif = ii - m.x1s;
                    // calculate density. If density < 0 or stopping time < dt: apply terminal velocity approximation.
                    m.dcons(specIND, IDN, kk, jj, ii) -= (dt / m.vol(kk, jj, ii) * (m.fdcons(specIND, IDN, 0, kkf, jjf, iif + 1) * m.f1a(kk, jj, ii + 1) - m.fdcons(specIND, IDN, 0, kkf, jjf, iif) * m.f1a(kk, jj, ii))
                                                        + dt / m.vol(kk, jj, ii) * (m.fdcons(specIND, IDN, 1, kkf, jjf + 1, iif) * m.f2a(kk, jj + 1, ii) - m.fdcons(specIND, IDN, 1, kkf, jjf, iif) * m.f2a(kk, jj, ii))
                                                        + dt / m.vol(kk, jj, ii) * (m.fdcons(specIND, IDN, 2, kkf + 1, jjf, iif) * m.f3a(kk + 1, jj, ii) - m.fdcons(specIND, IDN, 2, kkf, jjf, iif) * m.f3a(kk, jj, ii)));
                    bool denl0 = m.dcons(specIND, IDN, kk, jj, ii) < m.dminDensity;
                    bool stldt = m.stoppingtimemesh(specIND, kk, jj, ii) < dt;
                    if (denl0 || stldt){
                        if (denl0){
                            m.dcons(specIND, IDN, kk, jj, ii) = m.dminDensity;
                        }
                        double rhogradphix1;
                        double rhogradphix2;
                        double rhogradphix3;
                        #ifdef ENABLE_GRAVITY
                        rhogradphix1 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x1surface(kk, jj, ii + 1) - m.grav->Phi_grav_x1surface(kk, jj, ii)) / m.dx1p(kk, jj, ii);
                        rhogradphix2 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x2surface(kk, jj + 1, ii) - m.grav->Phi_grav_x2surface(kk, jj, ii)) / m.dx2p(kk, jj, ii);
                        rhogradphix3 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x3surface(kk + 1, jj, ii) - m.grav->Phi_grav_x3surface(kk, jj, ii)) / m.dx3p(kk, jj, ii);
                        #else   // set gravity to zero
                        rhogradphix1 = 0;
                        rhogradphix2 = 0;
                        rhogradphix3 = 0;
                        #endif // ENABLE_GRAVITY
                        dust_terminalvelocityapprixmation_rtp(m.prim(IV1, kk, jj, ii), m.prim(IV2, kk, jj, ii), m.prim(IV3, kk, jj, ii),
                                                              rhogradphix1, rhogradphix2, rhogradphix3,
                                                              m.dcons(specIND, IDN, kk, jj, ii), m.stoppingtimemesh(specIND, kk, jj, ii), m.x1v(ii), m.geo_cot(jj),
                                                              m.dcons(specIND, IM1, kk, jj, ii), m.dcons(specIND, IM2, kk, jj, ii), m.dcons(specIND, IM3, kk, jj, ii)
                                                              );
                        #ifdef DEBUG
                        if (denl0){
                            std::cout << m.dcons(specIND, IM1, kk, jj, ii) << '\t' << m.dcons(specIND, IM2, kk, jj, ii) << '\t' << m.dcons(specIND, IM3, kk, jj, ii) << std::endl << flush;
                        }
                        #endif // DEBUG
                    }
                    else{
                        // if density is fine, then calculate everything self-consistantly.
                        for (int dconsIND = 1; dconsIND < NUMCONS - 1; dconsIND++){
                            m.dcons(specIND, dconsIND, kk, jj, ii) -= (dt / m.vol(kk, jj, ii) * (m.fdcons(specIND, dconsIND, 0, kkf, jjf, iif + 1) * m.f1a(kk, jj, ii + 1) - m.fdcons(specIND, dconsIND, 0, kkf, jjf, iif) * m.f1a(kk, jj, ii))
                                                                     + dt / m.vol(kk, jj, ii) * (m.fdcons(specIND, dconsIND, 1, kkf, jjf + 1, iif) * m.f2a(kk, jj + 1, ii) - m.fdcons(specIND, dconsIND, 1, kkf, jjf, iif) * m.f2a(kk, jj, ii))
                                                                     + dt / m.vol(kk, jj, ii) * (m.fdcons(specIND, dconsIND, 2, kkf + 1, jjf, iif) * m.f3a(kk + 1, jj, ii) - m.fdcons(specIND, dconsIND, 2, kkf, jjf, iif) * m.f3a(kk, jj, ii)));
                        }
                        // geometry term

                        m.dcons(specIND, IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * \
                                    (m.dprim(specIND, IDN, kk, jj, ii) * (pow(m.dprim(specIND, IV2, kk, jj, ii), 2) + pow(m.dprim(specIND, IV3, kk, jj, ii), 2)));
                        m.dcons(specIND, IM2, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.dprim(specIND, IDN, kk, jj, ii) * m.dprim(specIND, IV1, kk, jj, ii) * m.dprim(specIND, IV2, kk, jj, ii);
                        m.dcons(specIND, IM2, kk, jj, ii) += dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.dprim(specIND, IDN, kk, jj, ii) * pow(m.dprim(specIND, IV3, kk, jj, ii), 2));
                        m.dcons(specIND, IM3, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.dprim(specIND, IDN, kk, jj, ii) * m.dprim(specIND, IV1, kk, jj, ii) * m.dprim(specIND, IV3, kk, jj, ii);
                        m.dcons(specIND, IM3, kk, jj, ii) -= dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.dprim(specIND, IDN, kk, jj, ii) * (m.dprim(specIND, IV2, kk, jj, ii) * m.dprim(specIND, IV3, kk, jj, ii)));

                        // source terms
                        double rhogradphix1;
                        double rhogradphix2;
                        double rhogradphix3;
                        #ifdef ENABLE_GRAVITY
                        rhogradphix1 = m.dprim(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x1surface(kk, jj, ii + 1) - m.grav->Phi_grav_x1surface(kk, jj, ii)) / m.dx1p(kk, jj, ii);
                        rhogradphix2 = m.dprim(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x2surface(kk, jj + 1, ii) - m.grav->Phi_grav_x2surface(kk, jj, ii)) / m.dx2p(kk, jj, ii);
                        rhogradphix3 = m.dprim(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x3surface(kk + 1, jj, ii) - m.grav->Phi_grav_x3surface(kk, jj, ii)) / m.dx3p(kk, jj, ii);
                        #else   // set gravity to zero
                        rhogradphix1 = 0;
                        rhogradphix2 = 0;
                        rhogradphix3 = 0;
                        #endif // ENABLE_GRAVITY

                        // apply gravity
                        m.dcons(specIND, IM1, kk, jj, ii) += rhogradphix1 * dt;
                        m.dcons(specIND, IM2, kk, jj, ii) += rhogradphix2 * dt;
                        m.dcons(specIND, IM3, kk, jj, ii) += rhogradphix3 * dt;
                        // gas drag
                        double vdust1 = m.dcons(specIND, IV1, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                        double vgas1  = m.prim(IV1, kk, jj, ii);
                        double vdust2 = m.dcons(specIND, IV2, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                        double vgas2  = m.prim(IV2, kk, jj, ii);
                        double vdust3 = m.dcons(specIND, IV3, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                        double vgas3  = m.prim(IV3, kk, jj, ii);
                        double rhodt_stime = m.dcons(specIND, IDN, kk, jj, ii) * dt / m.stoppingtimemesh(specIND, kk, jj, ii);
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
}
#endif


#if defined(CARTESIAN_COORD)
void advect_cons_dust_cartesian(mesh &m, double dt){
    // TODO: may need update
    #ifdef GPU
    #pragma acc parallel loop collapse (4) default(present)
    #else
    #pragma omp parallel for collapse (4) schedule (static)
    #endif
    for (int specIND  = 0; specIND < m.NUMSPECIES; specIND++){
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    int kkf = kk - m.x3s;
                    int jjf = jj - m.x2s;
                    int iif = ii - m.x1s;
                    m.dcons(specIND, IDN, kk, jj, ii)
                                -= (dt / m.dx1(ii) * (m.fdcons(specIND, IDN, 0, kkf, jjf, iif + 1) - m.fdcons(specIND, IDN, 0, kkf, jjf, iif))
                                  + dt / m.dx2(jj) * (m.fdcons(specIND, IDN, 1, kkf, jjf + 1, iif) - m.fdcons(specIND, IDN, 1, kkf, jjf, iif))
                                  + dt / m.dx3(kk) * (m.fdcons(specIND, IDN, 2, kkf + 1, jjf, iif) - m.fdcons(specIND, IDN, 2, kkf, jjf, iif)));
                    bool denl0 = m.dcons(specIND, IDN, kk, jj, ii) < m.dminDensity;
                    bool stldt = m.stoppingtimemesh(specIND, kk, jj, ii) < dt;
                    if (denl0 || stldt){
                        if (denl0){
                            m.dcons(specIND, IDN, kk, jj, ii) = m.dminDensity;
                        }
                        double rhogradphix1;
                        double rhogradphix2;
                        double rhogradphix3;
                        #ifdef ENABLE_GRAVITY
                        rhogradphix1 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x1surface(kk, jj, ii + 1) - m.grav->Phi_grav_x1surface(kk, jj, ii)) / m.dx1p(kk, jj, ii);
                        rhogradphix2 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x2surface(kk, jj + 1, ii) - m.grav->Phi_grav_x2surface(kk, jj, ii)) / m.dx2p(kk, jj, ii);
                        rhogradphix3 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x3surface(kk + 1, jj, ii) - m.grav->Phi_grav_x3surface(kk, jj, ii)) / m.dx3p(kk, jj, ii);
                        #else   // set gravity to zero
                        rhogradphix1 = 0;
                        rhogradphix2 = 0;
                        rhogradphix3 = 0;
                        #endif // ENABLE_GRAVITY
                        dust_terminalvelocityapprixmation_xyz(m.prim(IV1, kk, jj, ii), m.prim(IV2, kk, jj, ii), m.prim(IV3, kk, jj, ii),
                                                              rhogradphix1, rhogradphix2, rhogradphix3,
                                                              m.dcons(specIND, IDN, kk, jj, ii), m.stoppingtimemesh(specIND, kk, jj, ii),
                                                              m.dcons(specIND, IM1, kk, jj, ii), m.dcons(specIND, IM2, kk, jj, ii), m.dcons(specIND, IM3, kk, jj, ii)
                                                              );
                        #ifdef DEBUG
                        if (denl0){
                            std::cout << m.dcons(specIND, IM1, kk, jj, ii) << '\t' << m.dcons(specIND, IM2, kk, jj, ii) << '\t' << m.dcons(specIND, IM3, kk, jj, ii) << std::endl << flush;
                        }
                        #endif // DEBUG
                    }
                    else{
                        // if density is fine, then calculate everything self-consistantly.
                        for (int dconsIND = 1; dconsIND < NUMCONS - 1; dconsIND++){
                            m.dcons(specIND, dconsIND, kk, jj, ii) -= (dt / m.dx1(ii) * (m.fdcons(specIND, dconsIND, 0, kkf, jjf, iif + 1) - m.fdcons(specIND, dconsIND, 0, kkf, jjf, iif))
                                                                     + dt / m.dx2(jj) * (m.fdcons(specIND, dconsIND, 1, kkf, jjf + 1, iif) - m.fdcons(specIND, dconsIND, 1, kkf, jjf, iif))
                                                                     + dt / m.dx3(kk) * (m.fdcons(specIND, dconsIND, 2, kkf + 1, jjf, iif) - m.fdcons(specIND, dconsIND, 2, kkf, jjf, iif)));
                        }
                        // geometry term

                        // source terms
                        double rhogradphix1;
                        double rhogradphix2;
                        double rhogradphix3;
                        #ifdef ENABLE_GRAVITY
                        rhogradphix1 = m.dprim(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x1surface(kk, jj, ii + 1) - m.grav->Phi_grav_x1surface(kk, jj, ii)) / m.dx1p(kk, jj, ii);
                        rhogradphix2 = m.dprim(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x2surface(kk, jj + 1, ii) - m.grav->Phi_grav_x2surface(kk, jj, ii)) / m.dx2p(kk, jj, ii);
                        rhogradphix3 = m.dprim(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x3surface(kk + 1, jj, ii) - m.grav->Phi_grav_x3surface(kk, jj, ii)) / m.dx3p(kk, jj, ii);
                        #else   // set gravity to zero
                        rhogradphix1 = 0;
                        rhogradphix2 = 0;
                        rhogradphix3 = 0;
                        #endif // ENABLE_GRAVITY

                        // apply gravity
                        m.dcons(specIND, IM1, kk, jj, ii) += rhogradphix1 * dt;
                        m.dcons(specIND, IM2, kk, jj, ii) += rhogradphix2 * dt;
                        m.dcons(specIND, IM3, kk, jj, ii) += rhogradphix3 * dt;
                        // gas drag
                        double vdust1 = m.dcons(specIND, IV1, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                        double vgas1  = m.prim(IV1, kk, jj, ii);
                        double vdust2 = m.dcons(specIND, IV2, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                        double vgas2  = m.prim(IV2, kk, jj, ii);
                        double vdust3 = m.dcons(specIND, IV3, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                        double vgas3  = m.prim(IV3, kk, jj, ii);
                        double rhodt_stime = m.dcons(specIND, IDN, kk, jj, ii) * dt / m.stoppingtimemesh(specIND, kk, jj, ii);
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
}
#endif

void advect_cons_dust(mesh &m, double dt){
    #if defined(CARTESIAN_COORD)
    advect_cons_dust_cartesian(m, dt);
    #elif defined(SPHERICAL_POLAR_COORD)
    advect_cons_dust_sphericalpolar(m, dt);
    #else
        # error need coordinate defined
    #endif
}
