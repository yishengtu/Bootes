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
    double Eint;

    void setup(mesh &m, input_file &finput){
        init_unifdensity = 1.68322e-5; // 1e-14;         // uniform density initially
        central_point_mass = 1; // 1.989e33;      // Similar to the uniform case, this can be read from the input file.
        Eint = 1.0;
        //double temp_slope = (kT_mu_up - kT_mu_low) / (m.x1l - m.x1s);
        //double temp_const = kT_mu_low - m.x3s * temp_slope;
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

        for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
            for (int kk = m.x3s; kk < m.x3l; kk++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    for (int ii = m.x1s; ii < m.x1l; ii++){
                        double r = m.x1v(ii);
                        double OmegaK = sqrt((m.pconst.G * central_point_mass) / pow(r, 3));
                        double rho = init_unifdensity;
                        
                        m.dcons(specIND, IDN, kk, jj, ii) = m.cons(IDN, kk, jj, ii)/100.; // (double) specIND + 1.;
                        m.dcons(specIND, IM1, kk, jj, ii) = 0.0;
                        m.dcons(specIND, IM2, kk, jj, ii) = m.cons(IM2, kk, jj, ii)/100.;
                        m.dcons(specIND, IM3, kk, jj, ii) = 0.0; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                    }
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
                    double z = m.x2v(ii);
                    m.grav->grav_x1(kk, jj, ii) = - m.pconst.G * central_point_mass * r / pow(r * r + z * z, 1.5);
                    m.grav->grav_x2(kk, jj, ii) = 0;
                    m.grav->grav_x3(kk, jj, ii) = - m.pconst.G * central_point_mass * z / pow(r * r + z * z, 1.5);
                }
            }
        }

    }

    void setup_dust(mesh &m, input_file &finput){
        double smin    = finput.getDouble("smin");
        double smax    = finput.getDouble("smax");
        double rhodm   = finput.getDouble("srho");
        double length_scale = finput.getDouble("length_scale");
        int ns         = finput.getInt("num_species");
        m.NUMSPECIES = ns;
        m.rhodm = rhodm;

        m.GrainEdgeList = logspace(log10(smin), log10(smax), ns + 1, true);
        m.GrainSizeList.NewBootesArray(ns);
        for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
            double s2 = m.GrainEdgeList(specIND + 1);
            double s1 = m.GrainEdgeList(specIND);
            m.GrainSizeList(specIND) = pow((pow(s2, 4) - pow(s1, 4)) / (4 * (s2 - s1)), 1./3.);
            m.GrainSizeList(specIND) = m.GrainSizeList(specIND)/length_scale;
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
