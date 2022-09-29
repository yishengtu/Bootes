#include <cstdlib>
#include <cmath>

#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include "../algorithm/inoutput/input.hpp"
#include "../algorithm/gravity/gravity.hpp"
#include "../algorithm/index_def.hpp"


/** Setup a shearing box in 3D Cartesian Coordinate **/
namespace {
    double init_unifdensity;
    double kT_mu_up;
    double kT_mu_low;
    double central_point_mass = 0.3;
    double amin;      // m.GrainEdgeList(0);
    double amax;      // m.GrainEdgeList(m.NUMSPECIES);
    BootesArray<double> upperboundaryinitdustdensity;

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
        //double kT_mu = finput.getDouble("kT_mu");
        kT_mu_up = finput.getDouble("kT_mu_up");
        kT_mu_low = finput.getDouble("kT_mu_low");
        init_unifdensity = 1.6832940854700853e-4;
        double temp_slope = (kT_mu_up - kT_mu_low) / (m.x1l - m.x1s);
        double temp_const = kT_mu_low - m.x3s * temp_slope;
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double rho = init_unifdensity;
                    m.cons(IDN, kk, jj, ii) = rho;
                    m.cons(IM1, kk, jj, ii) = 0.0;    // 0.5 * sin(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0.0;    // 0.5 * cos(m.x2v(jj)) * m.cons(IDN, kk, jj, ii);
                    m.cons(IM3, kk, jj, ii) = 0.0;    // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                    // double kT_mu = kT_mu_up;
                    // m.cons(IEN, kk, jj, ii) = m.cons(IDN, kk, jj, ii) * temp / (m.hydro_gamma - 1.) / pow(m.x1v(ii), 0.5);
                    double temp = kT_mu_up; // + (1. - m.x3v(kk)) * 15736334.4567 ;
                    double IE = m.cons(IDN, kk, jj, ii) * temp / (m.hydro_gamma - 1.) / pow(m.x1v(ii), 6. / 7.);
                    m.cons(IEN, kk, jj, ii) = IE;
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
                                                m.hydro_gamma);
                }
            }
        }

        /** gravity **/
        double zero = 0.;
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
        // cout << "Alternatively, we can directly use the force grav_x1 etc." << endl << flush;

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

                    m.grav->grav_x1(kk, jj, ii) = - 0.09 * m.pconst.G * central_point_mass / (x * x); // + centfugal;
                    m.grav->grav_x2(kk, jj, ii) = 0;
                    m.grav->grav_x3(kk, jj, ii) = - OmegaKsq * z;        // Armitage equ. 234
                }
            }
        }

        // Right now, gravity is defined in main.cpp and time_integration.cpp.
        /** protection **/
        m.minTemp = kT_mu_low;
        // m.minDensity = 1e-4;

        // Dust
        #ifdef ENABLE_DUSTFLUID
        amin = finput.getDouble("ainimin");      // m.GrainEdgeList(0);
        amax = finput.getDouble("ainimax");      // m.GrainEdgeList(m.NUMSPECIES);
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
                            double rho = Adust * m.cons(IDN, kk, jj, ii) * (sqrt(amax) - sqrt(amin)) / (0.5 * pow(m.GrainSizeList(ss), 3)) * m.GrainMassList(ss);
                            m.dcons(ss, IDN, kk, jj, ii) = rho;
                            m.dcons(ss, IM1, kk, jj, ii) = rho * m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                            m.dcons(ss, IM2, kk, jj, ii) = rho * m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                            m.dcons(ss, IM3, kk, jj, ii) = rho * m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                        }
                    }
                }
            }
        }
        // Record the initial upper boundary dust density
        upperboundaryinitdustdensity.NewBootesArray(m.NUMSPECIES, m.nx2 + 2 * m.ng2, m.nx1 + 2 * m.ng1);
        for (int ss = 0; ss < m.NUMSPECIES; ss++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    upperboundaryinitdustdensity(ss, jj, ii) = m.dcons(ss, IDN, m.x3l - 1, jj, ii);
                }
            }
        }
        #endif // ENABLE_DUSTFLUID
    }


    void work_after_loop(mesh &m, double dt){
        // calculate accretion rate
        ;
        // re-calculate effective gravity based on situation now
        m.grav->zero_gravity(m);
        #ifdef GPU
        #pragma acc parallel loop collapse (3)
        #else
        #pragma omp parallel for collapse(3) schedule (static)
        #endif
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    // TODO: may not need the if statement
                    double x = m.x1v(ii);
                    double y = m.x2v(jj);
                    double z = m.x3v(kk);
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
                    double dpdr = (m.prim(IPN, kk, jj, indup) - m.prim(IPN, kk, jj, inddown)) / (m.x1v(indup) - m.x1v(inddown));
                    double eta = - 0.5 / m.cons(IDN, kk, jj, ii) * dpdr / gr;
                    double vphi = eta * sqrt(m.pconst.G * central_point_mass / x);
                    double centfugal = vphi * vphi / x;

                    m.grav->grav_x1(kk, jj, ii) = - 0.09 * m.pconst.G * central_point_mass / (x * x); // + centfugal;
                    m.grav->grav_x2(kk, jj, ii) = 0;

                    if (x > 0.6){
                        x = 0.6;
                    }
                    double OmegaKsq = (m.pconst.G * central_point_mass) / pow(x, 3);
                    m.grav->grav_x3(kk, jj, ii) = - OmegaKsq * z;        // Armitage equ. 234
                }
            }
        }

        // Put back in the pseudo-temperature profile
        #ifdef GPU
        #pragma acc parallel loop collapse (3)
        #else
        #pragma omp parallel for collapse(3) schedule (static)
        #endif
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    double temp = kT_mu_up; // + (1. - m.x3v(kk)) * 15736334.4567 ;
                    double KE = 0.5 * (pow(m.cons(IM1, kk, jj, ii), 2) + pow(m.cons(IM2, kk, jj, ii), 2) + pow(m.cons(IM3, kk, jj, ii), 2)) / m.cons(IDN, kk, jj, ii);
                    double IE = m.cons(IDN, kk, jj, ii) * temp / (m.hydro_gamma - 1.) / pow(m.x1v(ii), 6. / 7.);
                    // m.cons(IEN, kk, jj, ii) = std::max(m.cons(IEN, kk, jj, ii), KE + IE);
                    m.cons(IEN, kk, jj, ii) = KE + IE;
                }
            }
        }

        // Damp up-going z-direction sped
        #ifdef GPU
        #pragma acc parallel loop collapse (3)
        #else
        #pragma omp parallel for collapse(3) schedule (static)
        #endif
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    if (m.cons(IM3, kk, jj, ii) > 0){
                        m.cons(IM3, kk, jj, ii) *= 0.5;
                    }
                }
            }
        }
    }

    void apply_user_extra_boundary_condition(mesh &m){
        // upper boundary has fixed density of initial density
        #ifdef GPU
        #pragma acc parallel loop collapse (3)
        #else
        #pragma omp parallel for collapse(3) schedule (static)
        #endif
        for (int gind3 = 0; gind3 < m.ng3; gind3 ++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double v1 = m.prim(IV1, m.x3l + gind3, jj, ii);
                    double v2 = m.prim(IV2, m.x3l + gind3, jj, ii);
                    double v3 = m.prim(IV3, m.x3l + gind3, jj, ii);
                    double IE = m.cons(IEN, m.x3l + gind3, jj, ii) - 0.5 * m.cons(IDN, m.x3l + gind3, jj, ii) * (v1 * v1 + v2 * v2 + v3 * v3);
                    m.cons(IDN, m.x3l + gind3, jj, ii) = init_unifdensity;
                    m.cons(IM1, m.x3l + gind3, jj, ii) = v1 * m.cons(IDN, m.x3l + gind3, jj, ii);
                    m.cons(IM2, m.x3l + gind3, jj, ii) = v2 * m.cons(IDN, m.x3l + gind3, jj, ii);
                    m.cons(IM3, m.x3l + gind3, jj, ii) = v3 * m.cons(IDN, m.x3l + gind3, jj, ii);
                    m.cons(IEN, m.x3l + gind3, jj, ii) = IE + 0.5 * m.cons(IDN, m.x3l + gind3, jj, ii) * (v1 * v1 + v2 * v2 + v3 * v3);
                    m.prim(IDN, m.x3l + gind3, jj, ii) = m.cons(IDN, m.x3l + gind3, jj, ii);
                    m.prim(IPN, m.x3l + gind3, jj, ii) = pres(m.cons(IDN, m.x3l + gind3, jj, ii), m.cons(IEN, m.x3l + gind3, jj, ii),
                                                              m.cons(IM1, m.x3l + gind3, jj, ii), m.cons(IM2, m.x3l + gind3, jj, ii), m.cons(IM3, m.x3l + gind3, jj, ii), m.hydro_gamma);
                }
            }
        }

        // outer radial (x) boundary always has 0 radial speed.
        #ifdef GPU
        #pragma acc parallel loop collapse (3)
        #else
        #pragma omp parallel for collapse(3) schedule (static)
        #endif
        for (int gind1 = 0; gind1 < m.ng1; gind1 ++){
            for (int kk = m.x3s; kk < m.x3l; kk++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    //quan(IDN, kk, jj, x1l + gind1)     = quan(IDN, kk, jj, x1l - (gind1 + 1));
                    m.cons(IM1, kk, jj, m.x1l + gind1)     = 0.0;
                    m.prim(IV1, kk, jj, m.x1l + gind1)     = 0.0;
                    //quan(IM2, kk, jj, x1l + gind1)     = quan(IM2, kk, jj, x1l - (gind1 + 1));
                    //quan(IM3, kk, jj, x1l + gind1)     = quan(IM3, kk, jj, x1l - (gind1 + 1));
                    //quan(IEN, kk, jj, x1l + gind1)     = quan(IEN, kk, jj, x1l - (gind1 + 1));
                }
            }
        }

        // Dust
        #ifdef ENABLE_DUSTFLUID
        // upper boundary fixed density

        // upper boundary has fixed density of initial density
        #ifdef GPU
        #pragma acc parallel loop collapse (4)
        #else
        #pragma omp parallel for collapse(4) schedule (static)
        #endif
        for (int ss = 0; ss < m.NUMSPECIES; ss++){
            for (int gind3 = 0; gind3 < m.ng3; gind3 ++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    for (int ii = m.x1s; ii < m.x1l; ii++){
                        double v1 = m.dprim(ss, IV1, m.x3l + gind3, jj, ii);
                        double v2 = m.dprim(ss, IV2, m.x3l + gind3, jj, ii);
                        double v3 = m.dprim(ss, IV3, m.x3l + gind3, jj, ii);
                        m.dcons(ss, IDN, m.x3l + gind3, jj, ii) = upperboundaryinitdustdensity(ss, jj, ii);
                        m.dcons(ss, IM1, m.x3l + gind3, jj, ii) = v1 * m.dcons(ss, IDN, m.x3l + gind3, jj, ii);
                        m.dcons(ss, IM2, m.x3l + gind3, jj, ii) = v2 * m.dcons(ss, IDN, m.x3l + gind3, jj, ii);
                        m.dcons(ss, IM3, m.x3l + gind3, jj, ii) = v3 * m.dcons(ss, IDN, m.x3l + gind3, jj, ii);
                        m.dprim(ss, IDN, m.x3l + gind3, jj, ii) = m.dcons(IDN, m.x3l + gind3, jj, ii);
                    }
                }
            }
        }

        // TODO: outer radial boundary
        #ifdef GPU
        #pragma acc parallel loop collapse (4)
        #else
        #pragma omp parallel for collapse(4) schedule (static)
        #endif
        for (int ss = 0; ss < m.NUMSPECIES; ss++){
            for (int gind1 = 0; gind1 < m.ng1; gind1 ++){
                for (int kk = m.x3s; kk < m.x3l; kk++){
                    for (int jj = m.x2s; jj < m.x2l; jj++){
                        //quan(IDN, kk, jj, x1l + gind1)     = quan(IDN, kk, jj, x1l - (gind1 + 1));
                        m.dcons(ss, IM1, kk, jj, m.x1l + gind1)     = 0.0;
                        m.dprim(ss, IV1, kk, jj, m.x1l + gind1)     = 0.0;
                        //quan(IM2, kk, jj, x1l + gind1)     = quan(IM2, kk, jj, x1l - (gind1 + 1));
                        //quan(IM3, kk, jj, x1l + gind1)     = quan(IM3, kk, jj, x1l - (gind1 + 1));
                        //quan(IEN, kk, jj, x1l + gind1)     = quan(IEN, kk, jj, x1l - (gind1 + 1));
                    }
                }
            }
        }
        #endif // ENABLE_DUSTFLUID
    }


    void calculate_nu_vis(mesh &m){
        ;
    }
}
