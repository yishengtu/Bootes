#include "gravity.hpp"
#include "../BootesArray.hpp"
#include "../mesh/mesh.hpp"


gravity::gravity(){
    ;
}

void gravity::setup_Phimesh(int &tot_nx3, int &tot_nx2, int &tot_nx1){
    Phi_grav.NewBootesArray(tot_nx3, tot_nx2, tot_nx1);
    // TODO: can allocate a smaller array
    Phi_grav_x1surface.NewBootesArray(tot_nx3, tot_nx2, tot_nx1 + 1);
    Phi_grav_x2surface.NewBootesArray(tot_nx3, tot_nx2 + 1, tot_nx1);
    Phi_grav_x3surface.NewBootesArray(tot_nx3 + 1, tot_nx2, tot_nx1);
}


void gravity::pointsource_grav(mesh &m, double m_source, double x1_s, double x2_s, double x3_s){
    #pragma omp parallel for collapse (3)
    for (int kk = 0; kk < m.x3v.shape()[0]; kk ++){
        for (int jj = 0; jj < m.x2v.shape()[0]; jj ++){
            for (int ii = 0; ii < m.x1v.shape()[0]; ii ++){
                double x1_h = m.x1v(ii);
                double x2_h = m.x2v(jj);
                double x3_h = m.x3v(kk);
                double r;
                #if defined(CARTESIAN_COORD)
                    r = sqrt(pow(x1_h-x1_s, 2) + pow(x2_h-x2_s, 2) + pow(x3_h-x3_s, 2));
                #elif defined(SPHERICAL_POLAR_COORD)
                    r = sqrt(pow(x1_h, 2) + pow(x1_s, 2) - 2 * x1_h * x1_s * (sin(x2_h) * sin(x2_s) * cos(x3_h - x3_s) + cos(x2_h) * cos(x2_s)));
                #endif // defined
                Phi_grav(kk, jj, ii) = G * m_source / r;
            }
        }
    }
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
    if (m.dim > 1){
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l; kk ++){
            for (int jj = m.x2s; jj < m.x2l + 1; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    // TODO: may not need the if statement
                    if (ii != m.x1v.shape()[0] && kk != m.x3f.shape()[0]){
                        double fit_k = (Phi_grav(kk, jj, ii) - Phi_grav(kk, jj - 1, ii)) / (m.x2v(jj) - m.x2v(jj - 1));
                        double fit_b = Phi_grav(kk, jj, ii) - fit_k * m.x2v(jj);
                        Phi_grav_x2surface(kk, jj, ii) = fit_k * m.x2f(jj) + fit_b;
                    }
                }
            }
        }
    }
    else{
        Phi_grav_x2surface.set_uniform(0);
    }
    if (m.dim > 2){
        #pragma omp parallel for collapse (3)
        for (int kk = m.x3s; kk < m.x3l + 1; kk ++){
            for (int jj = m.x2s; jj < m.x2l; jj ++){
                for (int ii = m.x1s; ii < m.x1l; ii ++){
                    // TODO: may not need the if statement
                    if (ii != m.x1f.shape()[0] && jj != m.x2f.shape()[0]){
                        double fit_k = (Phi_grav(kk, jj, ii) - Phi_grav(kk - 1, jj, ii)) / (m.x3v(kk) - m.x3v(kk - 1));
                        double fit_b = Phi_grav(kk, jj, ii) - fit_k * m.x3v(kk);
                        Phi_grav_x3surface(kk, jj, ii) = fit_k * m.x3f(kk) + fit_b;
                    }
                }
            }
        }
    }
    else{
        Phi_grav_x3surface.set_uniform(0);
    }
}
