#ifndef MINMOD_RECONSTRUCT_HPP_
#define MINMOD_RECONSTRUCT_HPP_

#include "../BootesArray.hpp"
#include "../eos/momentum.hpp"
#include "../../defs.hpp"

void minmod(float &quanp1, float &quan, float &quanm1, float &dx_axis, float &dt, float &Vui, float &acs, float &BquanL, float &BquanR){
    float w = 0.0;
    float Dim = quan - quanm1; // Toro 13.28
    float Dip = quanp1 - quan;
    float Di = 0.5 * (1. + w) * Dim + 0.5 * (1. - w) * Dip;                   // Toro 13.27, 14.37
    float Dib; // 14.44
    float beta = 1.;
    if (Dip > 0.){
        Dib = max(max((float) 0., min(beta * Dim, Dip)), min(Dim, beta * Dip));
    }
    else{
        Dib = min(min((float) 0., max(beta * Dim, Dip)), max(Dim, beta * Dip));
    }
    float c = acs * dt / dx_axis;
    if (Vui > 0)     { ; }
    else if (Vui < 0){ c = -c; }
    else             { c = 0; }           // Vui = 0;
    BquanL = quan - 0.5 * (1 + c) * Dib;                // Toro 13.33
    BquanR = quan + 0.5 * (1 - c) * Dib;                // Toro 13.33
}


void reconstruct(mesh &m,
                   BootesArray<float> &valsL,
                   BootesArray<float> &valsR,
                   int &x1excess, int &x2excess, int &x3excess,
                   int &axis,
                   int &IMP,
                   float &dt
                   ){
    // Computation starts in first ghost zone, for first active cell left boundary flux
    float zero = 0;
    //cout << "---" << endl;
    #pragma omp parallel for collapse (3) schedule (static)
    for (int kk = -x3excess; kk < m.nx3 + x3excess; kk++){
        for (int jj = -x2excess; jj < m.nx2 + x2excess; jj++){
            for (int ii = -x1excess; ii < m.nx1 + x1excess; ii++){
                float dx_axis, a;
                float cs = soundspeed(m.cons(IDN, m.x3s + kk, m.x2s + jj, m.x1s + ii), m.prim(IPN, m.x3s + kk, m.x2s + jj, m.x1s + ii), m.hydro_gamma);
                #if defined(CARTESIAN_COORD)
                    if      (axis == 0) { dx_axis = m.dx1(ii+x1excess); a = max(cs + m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                    else if (axis == 1) { dx_axis = m.dx2(jj+x2excess); a = max(cs + m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                    else                { dx_axis = m.dx3(kk+x3excess); a = max(cs + m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                #elif defined(SPHERICAL_POLAR_COORD)
                    if      (axis == 0) {
                        dx_axis = m.dx1p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = max(cs + m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                    else if (axis == 1) {
                        dx_axis = m.dx2p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = max(cs + m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                    else                {
                        dx_axis = m.dx3p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = max(cs + m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                #else
                    # error need coordinate defined
                #endif
                //if (dx_axis == 0){ cout << axis << '\t' << kk << '\t' << jj << '\t' << ii << '\t' << dx_axis << endl << flush; }
                /** speeds **/
                float Vui   = vel(m.cons(IMP, m.x3s + kk, m.x2s + jj, m.x1s),
                                  m.cons(IDN, m.x3s + kk, m.x2s + jj, m.x1s)
                                  );
                /** conservatives **/
                //cout << " 1---" << endl << flush;
                float BrhoL, BrhoR;
                minmod(m.cons(IDN, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IDN, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                       m.cons(IDN, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       BrhoL,
                       BrhoR);
                float Bm1L, Bm1R;
                minmod(m.cons(IM1, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IM1, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                       m.cons(IM1, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       Bm1L,
                       Bm1R);
                float Bm2L, Bm2R;
                minmod(m.cons(IM2, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IM2, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                       m.cons(IM2, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       Bm2L,
                       Bm2R);
                float Bm3L, Bm3R;
                minmod(m.cons(IM3, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IM3, m.x3s + kk, m.x2s + jj, m.x1s + ii),
                       m.cons(IM3, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       Bm3L,
                       Bm3R);
                float BeneL, BeneR;
                minmod(m.cons(IEN, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IEN, m.x3s + kk, m.x2s + jj, m.x1s + ii),
                       m.cons(IEN, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       BeneL,
                       BeneR);

                // Left of a cell is the right of an edge.
                if (kk == -1 || jj == -1 || ii == -1){
                    ;
                }
                else {
                    valsR(axis, IDN, kk, jj, ii) = BrhoL;
                    valsR(axis, IM1, kk, jj, ii) = Bm1L;
                    valsR(axis, IM2, kk, jj, ii) = Bm2L;
                    valsR(axis, IM3, kk, jj, ii) = Bm3L;
                    valsR(axis, IEN, kk, jj, ii) = BeneL;
                }

                if (kk == m.nx3 || jj == m.nx2 || ii == m.nx1){
                    ;
                }
                else{
                    valsL(axis, IDN, kk + x3excess, jj + x2excess, ii + x1excess) = BrhoR;
                    valsL(axis, IM1, kk + x3excess, jj + x2excess, ii + x1excess) = Bm1R;
                    valsL(axis, IM2, kk + x3excess, jj + x2excess, ii + x1excess) = Bm2R;
                    valsL(axis, IM3, kk + x3excess, jj + x2excess, ii + x1excess) = Bm3R;
                    valsL(axis, IEN, kk + x3excess, jj + x2excess, ii + x1excess) = BeneR;
                }
            }
        }
    }
    //cout << "===" << endl << flush;
}


#endif
