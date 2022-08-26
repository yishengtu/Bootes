#include <cstdlib>
#include <cmath>

#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include "../algorithm/inoutput/input.hpp"
#include "../algorithm/gravity/gravity.hpp"
#include "../algorithm/index_def.hpp"


/** Setup a shearing box in 3D Cartesian Coordinate **/

void setup_dust(mesh &m, input_file &finput){
    double smin    = finput.getDouble("smin");
    double smax    = finput.getDouble("smax");
    double rhodm   = finput.getDouble("srho");
    int ns         = finput.getInt("num_species");
    m.NUMSPECIES = ns;
    m.rhodm = rhodm;

    m.GrainEdgeList = logspace(log10(smin), log10(smax), ns + 1, true);
    m.GrainSizeList.NewBootesArray(ns);
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
        double s2 = m.GrainEdgeList(specIND + 1);
        double s1 = m.GrainEdgeList(specIND);
        m.GrainSizeList(specIND) = pow((pow(s2, 4) - pow(s1, 4)) / (4 * (s2 - s1)), 1./3.);
    }
}


void setup(mesh &m, input_file &finput){
    double kT_mu = finput.getDouble("kT_mu");
    double asq = m.hydro_gamma * kT_mu;
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                double rho = 1.6832940854700853e-4;
                m.cons(IDN, kk, jj, ii) = rho;
                m.cons(IM1, kk, jj, ii) = 0.0; // 0.5 * sin(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                m.cons(IM2, kk, jj, ii) = 0.0; // 0.5 * cos(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                m.cons(IM3, kk, jj, ii) = 0.0; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                m.cons(IEN, kk, jj, ii) = rho * kT_mu / (m.hydro_gamma - 1.);
            }
        }
    }
    /** use Nakagawa et al. equ. 1.9 & 2.22 to calculate centrifugal force **/
    BootesArray<double> pressure; pressure.NewBootesArray(m.cons.shape()[1], m.cons.shape()[2], m.cons.shape()[3]);
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                pressure(kk, jj, ii) = pres(m.cons(IDN, kk, jj, ii), m.cons(IEN, kk, jj, ii),
                                            m.cons(IM1, kk, jj, ii), m.cons(IM2, kk, jj, ii), m.cons(IM3, kk, jj, ii),
                                            m.hydro_gamma);          // Armitage equ. 234
            }
        }
    }

    /** gravity **/
    double zero = 0.;
    double central_point_mass = 0.03;
    m.grav->zero_gravity(m);
    // z-direction gravity
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                double rho = 1.0;
                double x = m.x1v(ii);
                double y = m.x2v(jj);
                double z = m.x3v(kk);
                double OmegaKsq = (m.pconst.G * central_point_mass) / pow(x, 3);
                //m.grav->Phi_grav(kk, jj, ii) -= 0.5 * OmegaKsq * z * z;        // Armitage equ. 234
            }
        }
    }
    // x-direction centrifugal force
    cout << "Alternatively, we can directly use the force grav_x1 etc." << endl << flush;

    #pragma omp parallel for collapse (3)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                // TODO: may not need the if statement
                double x = m.x1v(ii);
                double y = m.x2v(jj);
                double z = m.x3v(kk);
                double OmegaKsq = (m.pconst.G * central_point_mass) / pow(x, 3);
                // Nakagawa et al. 1986
                double gr   = m.pconst.G * central_point_mass / (x * x);
                int indup;
                int inddown;
                if (ii == m.x1s){
                    inddown = ii;
                }
                else{
                    inddown = ii - 1;
                }
                if (ii == m.x1l - 1){
                    indup = ii;
                }
                else{
                    indup = ii + 1;
                }
                double dpdr = (pressure(kk, jj, indup) - pressure(kk, jj, inddown)) / (m.x1v(indup) - m.x1v(inddown));
                double eta = - 0.5 / m.cons(IDN, kk, jj, ii) * dpdr / gr;
                double vphi = eta * sqrt(m.pconst.G * central_point_mass / x);
                double centfugal = vphi * vphi / x;

                m.grav->grav_x1(kk, jj, ii) = - m.pconst.G * central_point_mass / (x * x) + centfugal;
                m.grav->grav_x2(kk, jj, ii) = 0;
                m.grav->grav_x3(kk, jj, ii) = - OmegaKsq * z;        // Armitage equ. 234
            }
        }
    }

    // Old: trying to construct gravity with potential
    /* BootesArray<double> centrifugal_potential;
    centrifugal_potential.NewBootesArray(m.cons.shape()[1], m.cons.shape()[2]);
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                double rho = 1.0;
                double x = m.x1v(ii);
                double y = m.x2v(jj);
                double z = m.x3v(kk);
                double OmegaKsq = (m.pconst.G * central_point_mass) / pow(x, 3);
                // Nakagawa et al. 1986
                double gr   = m.pconst.G * central_point_mass / (x * x);
                int indup;
                int inddown;
                if (ii == m.x1s){
                    inddown = ii;
                }
                else{
                    inddown = ii - 1;
                }
                if (ii == m.x1l - 1){
                    indup = ii;
                }
                else{
                    indup = ii + 1;
                }
                double dpdr = (pressure(kk, jj, indup) - pressure(kk, jj, inddown)) / (m.x1v(indup) - m.x1v(inddown));
                double eta = - 0.5 / m.cons(IDN, kk, jj, ii) * dpdr / gr;
                double vphi = eta * sqrt(m.pconst.G * central_point_mass / x);
                double centfugal = vphi * vphi / x;

                centrifugal_potential(kk, jj) += centfugal * (m.x1f(ii + 1) - m.x1f(ii));
                m.grav->Phi_grav(kk, jj, ii) += centrifugal_potential(kk, jj);
            }
        }
    }
    m.grav->add_pointsource_grav(m, central_point_mass, zero, zero, zero);
    m.grav->boundary_grav(m);
    m.grav->calc_surface_vals(m);
    */
    // Right now, gravity is defined in main.cpp and time_integration.cpp.
    /** protection **/
    m.minTemp = kT_mu;
    // m.minDensity = 1e-4;
}


void work_after_loop(mesh &m, double &dt){
    // calculate accretion rate
    ;
}

void calculate_nu_vis(mesh &m){
    ;
    /*
    m.nu_vis.NewBootesArray(m.x3v.shape()[0], m.x2v.shape()[0], m.x1v.shape()[0]);
    // Uniform viscosity for now
    // m.nu_vis.set_uniform(0.00000000001);

    for (int kk = 0; kk < m.x3v.shape()[0]; kk++){
        for (int jj = 0; jj < m.x2v.shape()[0]; jj++){
            for (int ii = 0; ii < m.x1v.shape()[0]; ii++){
                m.nu_vis(kk, jj, ii) = 0.001 * sqrt(m.pconst.G * m.UserScalers(0) / pow(m.x1v(ii), 3)) * pow(m.x1v(ii), 2);
            }
        }
    }
    */
}
