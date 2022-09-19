#include "const_recon.hpp"

#include "../BootesArray.hpp"
#include "../eos/momentum.hpp"
#include "../eos/eos.hpp"
#include "../../defs.hpp"
#include "../index_def.hpp"
#include "../mesh/mesh.hpp"

void const_recon(double &quanp1, double &quan, double &quanm1, double &dx_axis, double &dt, double &Vui, double &acs, double &BquanL, double &BquanR){
    BquanL = quan;
    BquanR = quan;
}


void reconstruct_const(mesh &m,
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
    for (int kk = -x3excess; kk < m.nx3 + x3excess; kk++){
        for (int jj = -x2excess; jj < m.nx2 + x2excess; jj++){
            for (int ii = -x1excess; ii < m.nx1 + x1excess; ii++){
                double dx_axis, a;
                double cs = soundspeed(m.cons(IDN, m.x3s + kk, m.x2s + jj, m.x1s + ii), m.prim(IPN, m.x3s + kk, m.x2s + jj, m.x1s + ii), m.hydro_gamma);
                #if defined(CARTESIAN_COORD)
                    if      (axis == 0) { dx_axis = m.dx1(ii+x1excess); a = std::max(cs + m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                    else if (axis == 1) { dx_axis = m.dx2(jj+x2excess); a = std::max(cs + m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                    else                { dx_axis = m.dx3(kk+x3excess); a = std::max(cs + m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii));}
                #elif defined(SPHERICAL_POLAR_COORD)
                    if      (axis == 0) {
                        dx_axis = m.dx1p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = std::max(cs + m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                    else if (axis == 1) {
                        dx_axis = m.dx2p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = std::max(cs + m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                    else                {
                        dx_axis = m.dx3p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = std::max(cs + m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                #elif defined(CYLINDRICAL_COORD)
                // still needs to be converted to cyl coords. 
                    if      (axis == 0) {
                        dx_axis = m.dx1p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = std::max(cs + m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV1, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                    else if (axis == 1) {
                        dx_axis = m.dx2p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = std::max(cs + m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV2, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                    else                {
                        dx_axis = m.dx3p(kk + x3excess, jj + x2excess, ii+x1excess);
                        a = std::max(cs + m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii), cs - m.prim(IV3, m.x3s + kk, m.x2s + jj, m.x1s + ii));
                    }
                #else
                    # error need coordinate defined
                #endif
                //if (dx_axis == 0){ cout << axis << '\t' << kk << '\t' << jj << '\t' << ii << '\t' << dx_axis << endl << flush; }
                /** speeds **/
                double Vui   = vel(m.cons(IMP, m.x3s + kk, m.x2s + jj, m.x1s + ii),
                                  m.cons(IDN, m.x3s + kk, m.x2s + jj, m.x1s + ii)
                                  );
                /** conservatives **/
                double BrhoL, BrhoR;
                const_recon(m.cons(IDN, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IDN, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                       m.cons(IDN, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       BrhoL,
                       BrhoR);
                double Bm1L, Bm1R;
                const_recon(m.cons(IM1, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IM1, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                       m.cons(IM1, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       Bm1L,
                       Bm1R);
                double Bm2L, Bm2R;
                const_recon(m.cons(IM2, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IM2, m.x3s + kk,            m.x2s + jj,            m.x1s + ii),
                       m.cons(IM2, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       Bm2L,
                       Bm2R);
                double Bm3L, Bm3R;
                const_recon(m.cons(IM3, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
                       m.cons(IM3, m.x3s + kk, m.x2s + jj, m.x1s + ii),
                       m.cons(IM3, m.x3s + kk - x3excess, m.x2s + jj - x2excess, m.x1s + ii - x1excess),
                       dx_axis, dt,
                       Vui, a,
                       Bm3L,
                       Bm3R);
                double BeneL, BeneR;
                const_recon(m.cons(IEN, m.x3s + kk + x3excess, m.x2s + jj + x2excess, m.x1s + ii + x1excess),
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
}

