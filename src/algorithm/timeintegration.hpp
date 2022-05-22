#ifndef TIME_INTEGRATION_HPP_
#define TIME_INTEGRATION_HPP_

#include "reconstruct/minmod.hpp"
//#include "reconstruct/const_recon.hpp"
#include "time_step/time_step.hpp"
#include "BootesArray.hpp"
#include "util.hpp"
#include "hydro/hllc.hpp"
#include "hydro/hll.hpp"
#include "hydro/donercell.hpp"
#include "boundary_condition/apply_bc.hpp"

void first_order(mesh &m, double &dt){

    /* commented out, first order doesn't need intermediate steps.
    BootesArray<double> rho_a; rho_a.NewBootesArray(m.rho.shape()[0], m.rho.shape()[1], m.rho.shape()[2]);
    BootesArray<double> mo1_a; mo1_a.NewBootesArray(m.rho.shape()[0], m.rho.shape()[1], m.rho.shape()[2]);
    BootesArray<double> mo2_a; mo2_a.NewBootesArray(m.rho.shape()[0], m.rho.shape()[1], m.rho.shape()[2]);
    BootesArray<double> mo3_a; mo3_a.NewBootesArray(m.rho.shape()[0], m.rho.shape()[1], m.rho.shape()[2]);
    BootesArray<double> ene_a; ene_a.NewBootesArray(m.rho.shape()[0], m.rho.shape()[1], m.rho.shape()[2]);     */

    // First order integration
    // (axis, z, y, x)
    // Step 1: calculate flux
    BootesArray<double> fcons;
    fcons.NewBootesArray(NUMCONS, 3, m.cons.shape()[1] + 1, m.cons.shape()[2] + 1, m.cons.shape()[3] + 1);

    // store the reconstructed value
    // index: (advecting direction, quantity, kk, jj, ii)
    BootesArray<double> valsL;
    BootesArray<double> valsR;
    #if defined(ENABLE_GRAVITY)
    valsL.NewBootesArray(3, NUMCONS + 1, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    valsR.NewBootesArray(3, NUMCONS + 1, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    #else
    valsL.NewBootesArray(3, NUMCONS, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    valsR.NewBootesArray(3, NUMCONS, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    #endif
    for (int axis = 0; axis < m.dim; axis ++){
        // step 1.1: reconstruct left/right values
        int x1excess, x2excess, x3excess;
        int IMP;        // Index of the velocity used in this axis
        if      (axis == 0){ x1excess = 1; x2excess = 0; x3excess = 0; IMP = IM1;}
        else if (axis == 1){ x1excess = 0; x2excess = 1; x3excess = 0; IMP = IM2;}
        else if (axis == 2){ x1excess = 0; x2excess = 0; x3excess = 1; IMP = IM3;}
        else { cout << "axis > 3!!!" << endl << flush; throw 1; }

        reconstruct(m, valsL, valsR, x1excess, x2excess, x3excess, axis, IMP, dt);
        // step 1.2: solve the Riemann problem. Use HLL for now, update conservative vars
        #pragma omp parallel for collapse (3) schedule (static)
        for (int kk = 0; kk < m.nx3 + x3excess; kk ++){
            for (int jj = 0; jj < m.nx2 + x2excess; jj ++){
                for (int ii = 0; ii < m.nx1 + x1excess; ii ++){
                    double valL[5];
                    double valR[5];
                    double fxs[5];
                    valL[IDN] = valsL(axis, IDN, kk, jj, ii); valR[IDN] = valsR(axis, IDN, kk, jj, ii);
                    valL[IM1] = valsL(axis, IM1, kk, jj, ii); valR[IM1] = valsR(axis, IM1, kk, jj, ii);
                    valL[IM2] = valsL(axis, IM2, kk, jj, ii); valR[IM2] = valsR(axis, IM2, kk, jj, ii);
                    valL[IM3] = valsL(axis, IM3, kk, jj, ii); valR[IM3] = valsR(axis, IM3, kk, jj, ii);
                    valL[IEN] = valsL(axis, IEN, kk, jj, ii); valR[IEN] = valsR(axis, IEN, kk, jj, ii);
                    # if defined(ENABLE_GRAVITY)
                        hll_grav(valL, valR, fxs,
                                 IMP,                   // the momentum term to add pressure and gravity
                                 valsL(axis, IGN, kk, jj, ii),      // phi
                                 valsR(axis, IGN, kk, jj, ii),
                                 m.hydro_gamma
                                 );
                    # else
                        hll( valL, valR, fxs,
                             IMP,                   // the momentum term to add pressure; shift by one index (since first index is density)
                             m.hydro_gamma
                             );
                    #endif
                    fcons(IDN, axis, kk, jj, ii) = fxs[IDN];
                    fcons(IM1, axis, kk, jj, ii) = fxs[IM1];
                    fcons(IM2, axis, kk, jj, ii) = fxs[IM2];
                    fcons(IM3, axis, kk, jj, ii) = fxs[IM3];
                    fcons(IEN, axis, kk, jj, ii) = fxs[IEN];
                }
            }
        }
    }
    // Need to set values in the fcons to zeros.
    #pragma omp parallel for collapse (4) schedule (static)
    for (int axis = m.dim; axis < 3; axis ++){
        for (int kk = 0; kk < fcons.shape()[2]; kk ++){
            for (int jj = 0; jj < fcons.shape()[3]; jj ++){
                for (int ii = 0; ii < fcons.shape()[4]; ii ++){
                    fcons(IDN, axis, kk, jj, ii) = 0.0;
                    fcons(IM1, axis, kk, jj, ii) = 0.0;
                    fcons(IM2, axis, kk, jj, ii) = 0.0;
                    fcons(IM3, axis, kk, jj, ii) = 0.0;
                    fcons(IEN, axis, kk, jj, ii) = 0.0;
                }
            }
        }
    }

    // step 2: time integrate to update CONSERVATIVE variables, solve Riemann Problem
    // First order for now
    #if defined(CARTESIAN_COORD)
        #pragma omp parallel for collapse (3) schedule (static)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    int kkf = kk - m.x3s;
                    int jjf = jj - m.x2s;
                    int iif = ii - m.x1s;
                    for (int consIND = 0; consIND < NUMCONS; consIND++){
                        m.cons(consIND, kk, jj, ii) -= (dt / m.dx1(ii) * (fcons(consIND, 0, kkf, jjf, iif + 1) - fcons(consIND, 0, kkf, jjf, iif))
                                                      + dt / m.dx2(jj) * (fcons(consIND, 1, kkf, jjf + 1, iif) - fcons(consIND, 1, kkf, jjf, iif))
                                                      + dt / m.dx3(kk) * (fcons(consIND, 2, kkf + 1, jjf, iif) - fcons(consIND, 2, kkf, jjf, iif)));
                    }
                }
            }
        }
    #elif defined(SPHERICAL_POLAR_COORD)
        #pragma omp parallel for collapse (3) schedule (static)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    int kkf = kk - m.x3s;
                    int jjf = jj - m.x2s;
                    int iif = ii - m.x1s;
                    // first take care of the divergent terms
                    for (int consIND = 0; consIND < NUMCONS; consIND++){
                        m.cons(consIND, kk, jj, ii) -= (dt / m.vol(kk, jj, ii) * (fcons(consIND, 0, kkf, jjf, iif + 1) * m.f1a(kk, jj, ii + 1) - fcons(consIND, 0, kkf, jjf, iif) * m.f1a(kk, jj, ii))
                                                      + dt / m.vol(kk, jj, ii) * (fcons(consIND, 1, kkf, jjf + 1, iif) * m.f2a(kk, jj + 1, ii) - fcons(consIND, 1, kkf, jjf, iif) * m.f2a(kk, jj, ii))
                                                      + dt / m.vol(kk, jj, ii) * (fcons(consIND, 2, kkf + 1, jjf, iif) * m.f3a(kk + 1, jj, ii) - fcons(consIND, 2, kkf, jjf, iif) * m.f3a(kk, jj, ii)));
                    }
                    // geometry term

                    /*   // lower order, but easier to understand..
                    m.cons(IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * \
                                (m.prim(IDN, kk, jj, ii) * (pow(m.prim(IV2, kk, jj, ii), 2) + pow(m.prim(IV3, kk, jj, ii), 2)) + 2 * m.prim(IPN, kk, jj, ii));

                    m.cons(IM2, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.prim(IDN, kk, jj, ii) * m.prim(IV1, kk, jj, ii) * m.prim(IV2, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) += dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * pow(m.prim(IV3, kk, jj, ii), 2) + m.prim(IPN, kk, jj, ii));

                    m.cons(IM3, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.prim(IDN, kk, jj, ii) * m.prim(IV1, kk, jj, ii) * m.prim(IV3, kk, jj, ii);
                    m.cons(IM3, kk, jj, ii) -= dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * (m.prim(IV2, kk, jj, ii) * m.prim(IV3, kk, jj, ii)));
                    */

                    double rp = m.x1f(ii + 1);
                    double rm = m.x1f(ii);
                    m.cons(IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * (pow(m.prim(IV2, kk, jj, ii), 2) + pow(m.prim(IV3, kk, jj, ii), 2)) + 2 * m.prim(IPN, kk, jj, ii));

                    m.cons(IM2, kk, jj, ii) -= dt * m.dx1(ii) / m.rV(ii) * (rm * rm * valsL(0, IM2, kkf, jjf, iif) + rp * rp * valsR(0, IM2, kkf, jjf, iif + 1));
                    m.cons(IM2, kk, jj, ii) += dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * pow(m.prim(IV3, kk, jj, ii), 2) + m.prim(IPN, kk, jj, ii));

                    m.cons(IM3, kk, jj, ii) -= dt * m.dx1(ii) / m.rV(ii) * (rm * rm * valsL(0, IM3, kkf, jjf, iif) + rp * rp * valsR(0, IM3, kkf, jjf, iif + 1));
                    m.cons(IM3, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.geo_cot(jj) / (m.geo_sm(jj) + m.geo_sp(jj)) *
                                                            (m.geo_sm(jj) * valsL(1, IM2, kkf, jjf, iif) * valsL(1, IM3, kkf, jjf, iif)
                                                            + m.geo_sp(jj) * valsR(1, IM2, kkf, jjf + 1, iif) * valsR(1, IM3, kkf, jjf + 1, iif));

                    #if defined(ENABLE_GRAVITY)
                    m.cons(IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * 2 * m.Phi_grav(kk, jj, ii) * m.prim(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) += dt * m.geo_cot(jj) * m.one_orgeo(ii) * m.Phi_grav(kk, jj, ii) * m.prim(IDN, kk, jj, ii);
                    #endif

                }
            }
        }
    #else
        # error need coordinate defined
    #endif

    // step 3: add source terms
    ;
    // protections
    #if defined (PROTECTION_PROTECTION)
    #pragma omp parallel for collapse (3)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                if (isnan(m.cons(IDN, kk, jj, ii)) || m.cons(IDN, kk, jj, ii) < 1e-8){
                    m.cons(IDN, kk, jj, ii) = 1e-8;
                }
                if (isnan(m.cons(IM1, kk, jj, ii))){
                    m.cons(IM1, kk, jj, ii) = 0;
                }
                if (isnan(m.cons(IM2, kk, jj, ii))){
                    m.cons(IM2, kk, jj, ii) = 0;
                }
                if (isnan(m.cons(IM3, kk, jj, ii))){
                    m.cons(IM3, kk, jj, ii) = 0;
                }
                if (isnan(m.cons(IEN, kk, jj, ii)) || m.cons(IEN, kk, jj, ii) < 1e-8){
                    m.cons(IEN, kk, jj, ii) = 1e-8;
                }
            }
        }
    }
    #endif // defined(PROTECTION_PROTECTION)

    // step 3: use E.O.S. and relations to get primitive variables.
    m.cons_to_prim();
    // step 4: apply boundary conditions
    apply_boundary_condition(m);

}

#endif // TIME_INTEGRATION_HPP_
