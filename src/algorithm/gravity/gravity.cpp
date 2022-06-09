#include "gravity.hpp"
#include "../BootesArray.hpp"
#include "../mesh/mesh.hpp"
#include "../physical_constants.hpp"
#include "../index_def.hpp"
#include <cmath>


gravity::gravity(){
    ;
}

void gravity::setup_Phimesh(int &tot_nx3, int &tot_nx2, int &tot_nx1){
    Phi_grav.NewBootesArray(tot_nx3, tot_nx2, tot_nx1);
    // TODO: can allocate a smaller array
    Phi_grav_x1surface.NewBootesArray(tot_nx3, tot_nx2, tot_nx1 + 1);
    Phi_grav_x2surface.NewBootesArray(tot_nx3, tot_nx2 + 1, tot_nx1);
    Phi_grav_x3surface.NewBootesArray(tot_nx3 + 1, tot_nx2, tot_nx1);
    grav_x1.NewBootesArray(tot_nx3, tot_nx2, tot_nx1);
    grav_x2.NewBootesArray(tot_nx3, tot_nx2, tot_nx1);
    grav_x3.NewBootesArray(tot_nx3, tot_nx2, tot_nx1);
}


void gravity::zero_gravity(mesh &m){
    #pragma omp parallel for collapse (3)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                Phi_grav(kk, jj, ii) = 0;
            }
        }
    }
}


void gravity::add_pointsource_grav(mesh &m, double &m_source, double &x1_s, double &x2_s, double &x3_s){
    #pragma omp parallel for collapse (3)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                double x1_h = m.x1v(ii);
                double x2_h = m.x2v(jj);
                double x3_h = m.x3v(kk);
                double r;
                #if defined(CARTESIAN_COORD)
                    r = sqrt(pow(x1_h-x1_s, 2) + pow(x2_h-x2_s, 2) + pow(x3_h-x3_s, 2));
                #elif defined(SPHERICAL_POLAR_COORD)
                    r = sqrt(pow(x1_h, 2) + pow(x1_s, 2) - 2 * x1_h * x1_s * (sin(x2_h) * sin(x2_s) * cos(x3_h - x3_s) + cos(x2_h) * cos(x2_s)));
                #endif // defined
                Phi_grav(kk, jj, ii) += m.pconst.G * m_source / r;
            }
        }
    }
}


void gravity::add_self_grav(mesh &m){
    BootesArray<double> mass_in_shell;
    mass_in_shell.NewBootesArray(m.cons.shape()[3]);
    double mass_in_last_shell = 0;
    for (int ii = m.x1s; ii < m.x1l; ii ++){
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                mass_in_last_shell += m.cons(IDN, kk, jj, ii) * m.vol(kk, jj, ii);
            }
        }
        mass_in_shell(ii) = mass_in_last_shell;
    }
    std::cout << std::endl << std::flush;
    #pragma omp parallel for
    for (int ii = m.x1s; ii < m.x1l; ii ++){
        double x1_h = m.x1v(ii);
        //double x2_h = m.x2v(jj);
        //double x3_h = m.x3v(kk);
        double r = x1_h;
        double Phi_shell = 0;
        for (int iis = ii; iis < m.x1l; iis ++){
            Phi_shell += m.pconst.G * mass_in_shell(iis) / pow(m.x1v(iis), 2) * m.dx1p(iis);
        }
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                Phi_grav(kk, jj, ii) += Phi_shell;
            }
        }
    }
}


void gravity::boundary_grav(mesh &m){
    // boundaries
    #pragma omp parallel
    {
        #pragma omp for collapse (3) nowait
        for (int gind1 = 0; gind1 < m.ng1; gind1 ++){
            for (int kk = m.x3s; kk < m.x3l; kk++){
                for (int jj = m.x2s; jj < m.x2l; jj++){
                    Phi_grav(kk, jj, m.x1s - 1 - gind1) = Phi_grav(kk, jj, m.x1s + gind1);
                    Phi_grav(kk, jj, m.x1l + gind1)     = Phi_grav(kk, jj, m.x1l - (gind1 + 1));
                }
            }
        }
        #pragma omp for collapse (3) nowait
        for (int gind2 = 0; gind2 < m.ng2; gind2 ++){
            for (int kk = m.x3s; kk < m.x3l; kk++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    Phi_grav(kk, m.x2s - 1 - gind2, ii) = Phi_grav(kk, m.x2s + gind2, ii);
                    Phi_grav(kk, m.x2l + gind2, ii)     = Phi_grav(kk, m.x2l - (gind2 + 1), ii);
                }
            }
        }
        #pragma omp for collapse (3) nowait
        for (int gind3 = 0; gind3 < m.ng3; gind3 ++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    Phi_grav(m.x3s - 1 - gind3, jj, ii) = Phi_grav(m.x3s + gind3, jj, ii);
                    Phi_grav(m.x3l + gind3, jj, ii)     = Phi_grav(m.x3l - (gind3 + 1), jj, ii);
                }
            }
        }
    }
    #pragma omp barrier
}


void gravity::calc_surface_vals(mesh &m){
    #pragma omp parallel for collapse (3)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l + 1; ii ++){
                // TODO: may not need the if statement
                if (jj != m.x2v.shape()[0] && kk != m.x3v.shape()[0]){
                    double fit_k = (Phi_grav(kk, jj, ii) - Phi_grav(kk, jj, ii - 1)) / (m.x1v(ii) - m.x1v(ii - 1));
                    double fit_b = Phi_grav(kk, jj, ii) - fit_k * m.x1v(ii);
                    Phi_grav_x1surface(kk, jj, ii) = fit_k * m.x1f(ii) + fit_b;
                }
            }
        }
    }
    #pragma omp parallel for collapse (3)
    for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                // TODO: may not need the if statement
                if (jj != m.x2v.shape()[0] && kk != m.x3v.shape()[0]){
                    grav_x1(kk, jj, ii) = (Phi_grav_x1surface(kk, jj, ii + 1) - Phi_grav_x1surface(kk, jj, ii)) / m.dx1p(kk, jj, ii);
                }
            }
        }
    }
    if (m.dim > 1){
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l + 1; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    // TODO: may not need the if statement
                    if (ii != m.x1v.shape()[0] && kk != m.x3v.shape()[0]){
                        double fit_k = (Phi_grav(kk, jj, ii) - Phi_grav(kk, jj - 1, ii)) / (m.x2v(jj) - m.x2v(jj - 1));
                        double fit_b = Phi_grav(kk, jj, ii) - fit_k * m.x2v(jj);
                        Phi_grav_x2surface(kk, jj, ii) = fit_k * m.x2f(jj) + fit_b;
                    }
                }
            }
        }
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    // TODO: may not need the if statement
                    if (jj != m.x2v.shape()[0] && kk != m.x3v.shape()[0]){
                        grav_x2(kk, jj, ii) = (Phi_grav_x2surface(kk, jj + 1, ii) - Phi_grav_x2surface(kk, jj, ii)) / m.dx2p(kk, jj, ii);
                    }
                }
            }
        }
    }
    else{
        Phi_grav_x2surface.set_uniform(0);
        grav_x2.set_uniform(0);
    }
    if (m.dim > 2){
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l + 1; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    // TODO: may not need the if statement
                    if (ii != m.x1v.shape()[0] && jj != m.x2v.shape()[0]){
                        double fit_k = (Phi_grav(kk, jj, ii) - Phi_grav(kk - 1, jj, ii)) / (m.x3v(kk) - m.x3v(kk - 1));
                        double fit_b = Phi_grav(kk, jj, ii) - fit_k * m.x3v(kk);
                        Phi_grav_x3surface(kk, jj, ii) = fit_k * m.x3f(kk) + fit_b;
                    }
                }
            }
        }

        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    // TODO: may not need the if statement
                    if (jj != m.x2v.shape()[0] && kk != m.x3v.shape()[0]){
                        grav_x3(kk, jj, ii) = (Phi_grav_x3surface(kk + 1, jj, ii) - Phi_grav_x3surface(kk, jj, ii)) / m.dx3p(kk, jj, ii);
                    }
                }
            }
        }
    }
    else{
        Phi_grav_x3surface.set_uniform(0);
        grav_x3.set_uniform(0);
    }
}





