#include <cstdlib>
#include <cmath>

#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include "../algorithm/inoutput/input.hpp"
#include "../algorithm/gravity/gravity.hpp"
#include "../algorithm/index_def.hpp"


/** Setup a Keplerian rotating initial condition **/
namespace {
    double init_unifdensity;
    double init_unifinternal;
    double central_point_mass;


    void setup(mesh &m, input_file &finput){
        init_unifdensity = 1;         // uniform density initially
        init_unifinternal = 1;        // uniform internal energy initially (uniform pressure)
        central_point_mass = 10;      // Similar to the uniform case, this can be read from the input file.

        double temp_slope = (kT_mu_up - kT_mu_low) / (m.x1l - m.x1s);
        double temp_const = kT_mu_low - m.x3s * temp_slope;
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double rho = init_unifdensity;
                    double r = m.x1v(ii);
                    double OmegaK = sqrt((m.pconst.G * central_point_mass) / pow(r, 3));
                    m.cons(IDN, kk, jj, ii) = rho;
                    m.cons(IM1, kk, jj, ii) = 0.0;                     // radial
                    m.cons(IM2, kk, jj, ii) = rho * OmegaK * r;        // azimuthal
                    m.cons(IM3, kk, jj, ii) = 0.0;                     // height
                    m.cons(IEN, kk, jj, ii) = 0.5 * rho * pow(OmegaK * r, 2) + Eint;   // = KE + IE
                }
            }
        }

        /** gravity **/
        double zero = 0.;
        m.grav->zero_gravity(m); // first initialize the values in this array
        // Then put in gravity
        // Technically, the best way to do this is to input gravitational potential, then let the program calculate the gravitational acceleration
        // However, since we are only doing a static state, it is easier to put in acceleration directly.
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    double r = m.x1v(ii);
                    m.grav->grav_x1(kk, jj, ii) = - m.pconst.G * central_point_mass / (r * r);
                    m.grav->grav_x2(kk, jj, ii) = 0;
                    m.grav->grav_x3(kk, jj, ii) = 0;
                }
            }
        }

    }


    void work_after_loop(mesh &m, double &dt){
        /** gravity **/
        // put back in gravity (may not necessary but is safe)
        double zero = 0.;
        m.grav->zero_gravity(m); // first initialize the values in this array
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    double r = m.x1v(ii);
                    m.grav->grav_x1(kk, jj, ii) = - m.pconst.G * central_point_mass / (r * r);
                    m.grav->grav_x2(kk, jj, ii) = 0;
                    m.grav->grav_x3(kk, jj, ii) = 0;
                }
            }
        }
    }

    void apply_user_extra_boundary_condition(mesh &m){
        ;
    }


    void calculate_nu_vis(mesh &m){
        ;
    }
}
