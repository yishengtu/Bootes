#include <cstdlib>
#include <cmath>

#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include "../algorithm/inoutput/input.hpp"
#include "../algorithm/gravity/gravity.hpp"
#include "../algorithm/index_def.hpp"


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
    double A = finput.getDouble("A");
    double omega = finput.getDouble("omega");
    double asq = m.hydro_gamma * kT_mu;
    double d_coef = asq * A / (4 * M_PI * m.pconst.G);
    double e_coef = asq * d_coef / (m.hydro_gamma * (m.hydro_gamma - 1.));
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                double r = m.x1v(ii);
                // double a0 = finput.getDouble("iso_sound_speed");
                double rm = m.x1f(ii);
                double rp = m.x1f(ii + 1);
                double invrsq = 3 * (rp - rm) / (rp * rp * rp - rm * rm * rm);
                double sintheta = - (cos(m.x2f(jj + 1)) - cos(m.x2f(jj))) / (m.x2f(jj + 1) - m.x2f(jj));
                double rho = d_coef * invrsq;
                double vphi = omega * r * sintheta;
                m.cons(IDN, kk, jj, ii) = rho;
                m.cons(IM1, kk, jj, ii) = 0.0; // 0.5 * sin(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                m.cons(IM2, kk, jj, ii) = 0.0; // 0.5 * cos(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                //m.cons(IM1, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                //m.cons(IM2, kk, jj, ii) += 0.01 * ((double) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                m.cons(IM3, kk, jj, ii) = rho * vphi; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                m.cons(IEN, kk, jj, ii) = e_coef * invrsq + 0.5 * rho * vphi * vphi;
            }
        }
    }

    #ifdef ENABLE_GRAVITY
    /** gravity **/
    double point_mass = d_coef * 4 * M_PI * m.x1f(m.x1s);

    m.UserScalers.NewBootesArray(1);    // 1 element: the point mass at the center
    m.UserScalers(0) = point_mass;

    double zero = 0.;
    m.grav->zero_gravity(m);
    m.grav->add_pointsource_grav(m, m.UserScalers(0), zero, zero, zero);
    m.grav->add_self_grav(m);
    m.grav->boundary_grav(m);
    m.grav->calc_surface_vals(m);
    // Right now, gravity is defined in main.cpp and time_integration.cpp.
    #endif // ENABLE_GRAVITY
    /** protection **/
    m.minTemp = kT_mu;
    m.minDensity = 1e-4;

    #ifdef ENABLE_DUSTFLUID
    double amin = finput.getDouble("ainimin");      // m.GrainEdgeList(0);
    double amax = finput.getDouble("ainimax");      // m.GrainEdgeList(m.NUMSPECIES);
    m.dminDensity = finput.getDouble("dminDensity");

    double Adust = 3. / (100. * 8 * M_PI * m.rhodm * (sqrt(amax) - sqrt(amin)));

    for (int ss = 0; ss < m.NUMSPECIES; ss++){
        if (m.GrainEdgeList(ss + 1) < amin || m.GrainEdgeList(ss) > amax){
            for (int kk = m.x3s; kk < m.x3l; kk++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    for (int ii = m.x1s; ii < m.x1l; ii++){
                        m.dcons(ss, IDN, kk, jj, ii) = m.dminDensity;
                        m.dcons(ss, IM1, kk, jj, ii) = m.dminDensity * m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        m.dcons(ss, IM2, kk, jj, ii) = m.dminDensity * m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        m.dcons(ss, IM3, kk, jj, ii) = m.dminDensity * m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    }
                }
            }
        }
        else{
            for (int kk = m.x3s; kk < m.x3l; kk++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    for (int ii = m.x1s; ii < m.x1l; ii++){
                        double upperlim = min(amax, m.GrainEdgeList(ss + 1));
                        double lowerlim = max(amin, m.GrainEdgeList(ss));
                        double rho = Adust * m.cons(IDN, kk, jj, ii) * (sqrt(upperlim) - sqrt(lowerlim)) / (0.5 * pow(m.GrainSizeList(ss), 3)) * m.GrainMassList(ss);
                        m.dcons(ss, IDN, kk, jj, ii) = rho;
                        m.dcons(ss, IM1, kk, jj, ii) = rho * m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        m.dcons(ss, IM2, kk, jj, ii) = rho * m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        m.dcons(ss, IM3, kk, jj, ii) = rho * m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    }
                }
            }
        }
    }

    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                cout << m.dcons(0, IDN, kk, jj, ii) << '\t';
            }
        }
        cout << endl << flush;
    }

    #endif // ENABLE_DUSTFLUID
}


void work_after_loop(mesh &m, double &dt){
    // calculate accretion rate
    ;
    double dmdt = 0;
    for (int jj = m.x2s; jj < m.x2l; jj++){
        dmdt -= m.cons(IDN, 0, jj, m.x1s) * std::min(m.prim(IV1, 0, jj, m.x1s), (double) 0) * m.f1a(0, jj, m.x1s);
    }
    m.UserScalers(0) += dmdt * dt;
    m.grav->zero_gravity(m);
    double zero = 0.;
    m.grav->add_pointsource_grav(m, m.UserScalers(0), zero, zero, zero);
    m.grav->add_self_grav(m);
    m.grav->boundary_grav(m);
    m.grav->calc_surface_vals(m);
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
