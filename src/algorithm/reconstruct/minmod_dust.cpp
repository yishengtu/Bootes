#ifndef MINMOD_RECONSTRUCT_DUST_HPP_
#define MINMOD_RECONSTRUCT_DUST_HPP_

#include "minmod_dust.hpp"
#include "../BootesArray.hpp"
#include "../eos/momentum.hpp"
#include "../../defs.hpp"
#include "../index_def.hpp"
#include "../mesh/mesh.hpp"

void minmod_dust(double &quanp1, double &quan, double &quanm1, double &dx_axis, double &dt, double &Vui, double &acs, double &BquanL, double &BquanR){
    double w = 0.0;
    double Dim = quan - quanm1; // Toro 13.28
    double Dip = quanp1 - quan;
    double Di = 0.5 * (1. + w) * Dim + 0.5 * (1. - w) * Dip;                   // Toro 13.27, 14.37
    double Dib; // 14.44
    double beta = 1.;
    if (Dip > 0.){
        Dib = std::max(std::max((double) 0., std::min(beta * Dim, Dip)), std::min(Dim, beta * Dip));
    }
    else{
        Dib = std::min(std::min((double) 0., std::max(beta * Dim, Dip)), std::max(Dim, beta * Dip));
    }
    double c = acs * dt / dx_axis;
    if (Vui > 0)     { ; }
    else if (Vui < 0){ c = -c; }
    else             { c = 0; }           // Vui = 0;
    BquanL = quan - 0.5 * (1 + c) * Dib;                // Toro 13.33
    BquanR = quan + 0.5 * (1 - c) * Dib;                // Toro 13.33
}


void reconstruct_dust(mesh &m,
                   BootesArray<double> &valsL,
                   BootesArray<double> &valsR,
                   int &x1excess, int &x2excess, int &x3excess,
                   int &axis,
                   int &IMP,
                   double &dt
                   ){
    // Computation starts in first ghost zone, for first active cell left boundary flux
    double zero = 0;
    #pragma omp parallel for collapse (3) schedule (static)
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
        for (int kk = -x3excess; kk < m.nx3 + x3excess; kk++){
            for (int jj = -x2excess; jj < m.nx2 + x2excess; jj++){
                for (int ii = -x1excess; ii < m.nx1 + x1excess; ii++){
                    double dx_axis, a;
                    #if defined(CARTESIAN_COORD)
                        if      (axis == 0) { dx_axis = m.dx1(ii+x1excess); a = std::max(m.dprim(specIND, IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii), - m.dprim(specIND, IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                        else if (axis == 1) { dx_axis = m.dx2(jj+x2excess); a = std::max(m.dprim(specIND, IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii), - m.dprim(specIND, IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                        else                { dx_axis = m.dx3(kk+x3excess); a = std::max(m.dprim(specIND, IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii), - m.dprim(specIND, IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                    #elif defined(SPHERICAL_POLAR_COORD)
                        if      (axis == 0) {
                            dx_axis = m.dx1p(kk + x3excess, jj + x2excess, ii+x1excess);
                            a = std::max(m.dprim(specIND, IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii), - m.dprim(specIND, IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                        }
                        else if (axis == 1) {
                            dx_axis = m.dx2p(kk + x3excess, jj + x2excess, ii+x1excess);
                            a = std::max(m.dprim(specIND, IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii), - m.dprim(specIND, IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                        }
                        else                {
                            dx_axis = m.dx3p(kk + x3excess, jj + x2excess, ii+x1excess);
                            a = std::max(m.dprim(specIND, IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii), - m.dprim(specIND, IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                        }
                    #else
                        # error need coordinate defined
                    #endif

                    /** speeds **/
                    double Vui   = vel(m.dcons(specIND, IMP, m.x3s + kk, m.x2s + jj, m.x1s),
                                       m.dcons(specIND, IDN, m.x3s + kk, m.x2s + jj, m.x1s)
                                       );
                    /** dconservatives **/
                    double BrhoL, BrhoR;
                    minmod_dust(m.dcons(specIND, IDN, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                                m.dcons(specIND, IDN, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                                m.dcons(specIND, IDN, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                                dx_axis, dt,
                                Vui, a,
                                BrhoL,
                                BrhoR);
                    double Bm1L, Bm1R;
                    minmod_dust(m.dcons(specIND, IM1, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                                m.dcons(specIND, IM1, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                                m.dcons(specIND, IM1, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                                dx_axis, dt,
                                Vui, a,
                                Bm1L,
                                Bm1R);
                    double Bm2L, Bm2R;
                    minmod_dust(m.dcons(specIND, IM2, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                                m.dcons(specIND, IM2, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                                m.dcons(specIND, IM2, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                                dx_axis, dt,
                                Vui, a,
                                Bm2L,
                                Bm2R);
                    double Bm3L, Bm3R;
                    minmod_dust(m.dcons(specIND, IM3, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                                m.dcons(specIND, IM3, m.x3s + kk, m.x2s + jj, m.x1s + ii),
                                m.dcons(specIND, IM3, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                                dx_axis, dt,
                                Vui, a,
                                Bm3L,
                                Bm3R);

                    // Left of a cell is the right of an edge.
                    if (kk == -1 || jj == -1 || ii == -1){
                        ;
                    }
                    else {
                        valsR(specIND, axis, IDN, kk, jj, ii) = BrhoL;
                        valsR(specIND, axis, IM1, kk, jj, ii) = Bm1L;
                        valsR(specIND, axis, IM2, kk, jj, ii) = Bm2L;
                        valsR(specIND, axis, IM3, kk, jj, ii) = Bm3L;
                    }

                    if (kk == m.nx3 || jj == m.nx2 || ii == m.nx1){
                        ;
                    }
                    else{
                        valsL(specIND, axis, IDN, kk + x3excess, jj + x2excess, ii + x1excess) = BrhoR;
                        valsL(specIND, axis, IM1, kk + x3excess, jj + x2excess, ii + x1excess) = Bm1R;
                        valsL(specIND, axis, IM2, kk + x3excess, jj + x2excess, ii + x1excess) = Bm2R;
                        valsL(specIND, axis, IM3, kk + x3excess, jj + x2excess, ii + x1excess) = Bm3R;
                    }
                }
            }
        }
    }
}


#endif
