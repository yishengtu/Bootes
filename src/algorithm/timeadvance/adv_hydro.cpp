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


void calc_flux(mesh &m, double &dt, BootesArray<double> &fcons, BootesArray<double> &valsL, BootesArray<double> &valsR){
    // store the reconstructed value
    // index: (advecting direction, quantity, kk, jj, ii)
    for (int axis = 0; axis < m.dim; axis ++){
        // step 1.1: reconstruct left/right values
        int x1excess, x2excess, x3excess;
        int IMP;        // Index of the velocity used in this axis
        if      (axis == 0){ x1excess = 1; x2excess = 0; x3excess = 0; IMP = IM1;}
        else if (axis == 1){ x1excess = 0; x2excess = 1; x3excess = 0; IMP = IM2;}
        else if (axis == 2){ x1excess = 0; x2excess = 0; x3excess = 1; IMP = IM3;}
        else { cout << "axis > 3!!!" << endl << flush; throw 1; }

        reconstruct_minmod(m, valsL, valsR, x1excess, x2excess, x3excess, axis, IMP, dt);
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
                    #ifdef ENABLE_TEMPERATURE_PROTECTION
                    valL[IEN] = energy_from_temperature_protection(valL[IDN], valL[IEN], valL[IM1], valL[IM2], valL[IM3], m.minTemp, m.hydro_gamma);
                    valR[IEN] = energy_from_temperature_protection(valR[IDN], valR[IEN], valR[IM1], valR[IM2], valR[IM3], m.minTemp, m.hydro_gamma);
                    #endif // ENABLE_TEMPERATURE_PROTECTION
                    hlle(valL, valR, fxs,
                         IMP,                   // the momentum term to add pressure; shift by one index (since first index is density)
                         m.hydro_gamma
                         );
                    fcons(IDN, axis, kk, jj, ii) = fxs[IDN];
                    fcons(IM1, axis, kk, jj, ii) = fxs[IM1];
                    fcons(IM2, axis, kk, jj, ii) = fxs[IM2];
                    fcons(IM3, axis, kk, jj, ii) = fxs[IM3];
                    fcons(IEN, axis, kk, jj, ii) = fxs[IEN];
                }
            }
        }
    }
    // Need to set unused values in the fdcons to zeros.
    // To do so the axis goes from "number of active axis" to 3
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
}


void advect_cons(mesh &m, double &dt, BootesArray<double> &fcons, BootesArray<double> &valsL, BootesArray<double> &valsR){
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
                    // lower order, but easier to understand..
                    m.cons(IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * \
                                (m.prim(IDN, kk, jj, ii) * (pow(m.prim(IV2, kk, jj, ii), 2) + pow(m.prim(IV3, kk, jj, ii), 2)) + 2 * m.prim(IPN, kk, jj, ii));
                    m.cons(IM2, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.prim(IDN, kk, jj, ii) * m.prim(IV1, kk, jj, ii) * m.prim(IV2, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) += dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * pow(m.prim(IV3, kk, jj, ii), 2) + m.prim(IPN, kk, jj, ii));
                    m.cons(IM3, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.prim(IDN, kk, jj, ii) * m.prim(IV1, kk, jj, ii) * m.prim(IV3, kk, jj, ii);
                    m.cons(IM3, kk, jj, ii) -= dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * (m.prim(IV2, kk, jj, ii) * m.prim(IV3, kk, jj, ii)));

                    /*
                    m.cons(IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * (pow(m.prim(IV2, kk, jj, ii), 2) + pow(m.prim(IV3, kk, jj, ii), 2)) + 2 * m.prim(IPN, kk, jj, ii));

                    m.cons(IM2, kk, jj, ii) -= dt * m.dx1(ii) / m.rV(ii) * (m.rsq(ii) * valsL(0, IM2, kkf, jjf, iif) + m.rsq(ii + 1) * valsR(0, IM2, kkf, jjf, iif + 1));
                    m.cons(IM2, kk, jj, ii) += dt * m.geo_cot(jj) * m.one_orgeo(ii) * (m.prim(IDN, kk, jj, ii) * pow(m.prim(IV3, kk, jj, ii), 2) + m.prim(IPN, kk, jj, ii));

                    m.cons(IM3, kk, jj, ii) -= dt * m.dx1(ii) / m.rV(ii) * (m.rsq(ii) * valsL(0, IM3, kkf, jjf, iif) + m.rsq(ii + 1) * valsR(0, IM3, kkf, jjf, iif + 1));
                    m.cons(IM3, kk, jj, ii) -= dt * m.one_orgeo(ii) * m.geo_cot(jj) / (m.geo_sm(jj) + m.geo_sp(jj)) *
                                                            (m.geo_sm(jj) * valsL(1, IM2, kkf, jjf, iif) * valsL(1, IM3, kkf, jjf, iif)
                                                            + m.geo_sp(jj) * valsR(1, IM2, kkf, jjf + 1, iif) * valsR(1, IM3, kkf, jjf + 1, iif));
                    */
                }
            }
        }
    #elif defined(CYLINDRICAL_COORD)
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
                    m.cons(IM1, kk, jj, ii) += dt * m.one_orgeo(ii) * \
                                ((m.prim(IDN, kk, jj, ii) * pow(m.prim(IV2, kk, jj, ii), 2)) + m.prim(IPN, kk, jj, ii));
                    m.cons(IM2, kk, jj, ii) -= dt * m.one_orgeo(ii) * \
                                (m.prim(IDN, kk, jj, ii) * (m.prim(IV1, kk, jj, ii) * m.prim(IV2, kk, jj, ii)));
                    m.cons(IM3, kk, jj, ii) -= 0.0;
                    //m.cons(IEN, kk, jj, ii) += dt * m.one_orgeo(ii) * \
                    //            (m.prim(IDN, kk, jj, ii) * \
                    //            (pow(m.prim(IV3, kk, jj, ii),3) + m.prim(IV1, kk, jj, ii) * pow(m.prim(IV3, kk, jj, ii),2)));
                }
            }
        }
    #else
        # error need coordinate defined
    #endif
}


#ifdef DENSITY_PROTECTION
void protection(mesh &m, double &minDensity){
    #pragma omp parallel for collapse (3)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                if (isnan(m.cons(IDN, kk, jj, ii)) || m.cons(IDN, kk, jj, ii) < minDensity){
                    m.cons(IDN, kk, jj, ii) = minDensity;
                }
                if (isnan(m.cons(IM1, kk, jj, ii))){
                    m.cons(IM1, kk, jj, ii) = minDensity * m.prim(IM1, kk, jj, ii);
                }
                if (isnan(m.cons(IM2, kk, jj, ii))){
                    m.cons(IM2, kk, jj, ii) = minDensity * m.prim(IM2, kk, jj, ii);
                }
                if (isnan(m.cons(IM3, kk, jj, ii))){
                    m.cons(IM3, kk, jj, ii) = minDensity * m.prim(IM3, kk, jj, ii);
                }
                if (isnan(m.cons(IEN, kk, jj, ii)) || m.cons(IEN, kk, jj, ii) < minDensity){
                    m.cons(IEN, kk, jj, ii) = minDensity;
                }
            }
        }
    }
}
#endif // DENSITY_PROTECTION

#ifdef ENABLE_TEMPERATURE_PROTECTION
void temperature_protection(mesh &m, double &minTemp){
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                m.cons(IEN, kk, jj, ii) = energy_from_temperature_protection(m.cons(IDN, kk, jj, ii), m.cons(IEN, kk, jj, ii),
                                                                             m.cons(IM1, kk, jj, ii), m.cons(IM2, kk, jj, ii), m.cons(IM3, kk, jj, ii), minTemp, m.hydro_gamma);
            }
        }
    }
}
#endif // ENABLE_TEMPERATURE_PROTECTION


