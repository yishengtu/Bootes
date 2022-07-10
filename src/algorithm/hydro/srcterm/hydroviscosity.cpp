#include "hydrograv.hpp"

#include "../../reconstruct/minmod.hpp"
//#include "reconstruct/const_recon.hpp"
#include "../../BootesArray.hpp"
#include "../../util/util.hpp"
#include "../hll.hpp"
#include "../../index_def.hpp"
#include "../../mesh/mesh.hpp"


void calcgradV(mesh &m, BootesArray<double> &gradV){
    gradV.set_uniform(0);
    // Step 1: calculate speeds on the boundaries
    BootesArray<double> v1bound1; v1bound1.NewBootesArray(gradV.shape()[0], gradV.shape()[1], gradV.shape()[2] + 1);
    //BootesArray<double> v2bound1; v2bound1.NewBootesArray(gradV.shape()[0], gradV.shape()[1], gradV.shape()[2] + 1);
    //BootesArray<double> v3bound1; v3bound1.NewBootesArray(gradV.shape()[0], gradV.shape()[1], gradV.shape()[2] + 1);
    //BootesArray<double> v1bound2; v1bound2.NewBootesArray(gradV.shape()[0], gradV.shape()[1] + 1, gradV.shape()[2]);
    BootesArray<double> v2bound2; v2bound2.NewBootesArray(gradV.shape()[0], gradV.shape()[1] + 1, gradV.shape()[2]);
    //BootesArray<double> v3bound2; v3bound2.NewBootesArray(gradV.shape()[0], gradV.shape()[1] + 1, gradV.shape()[2]);
    //BootesArray<double> v1bound3; v1bound3.NewBootesArray(gradV.shape()[0] + 1, gradV.shape()[1], gradV.shape()[2]);
    //BootesArray<double> v2bound3; v2bound3.NewBootesArray(gradV.shape()[0] + 1, gradV.shape()[1], gradV.shape()[2]);
    BootesArray<double> v3bound3; v3bound3.NewBootesArray(gradV.shape()[0] + 1, gradV.shape()[1], gradV.shape()[2]);

    // x1-boundaries
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l + 1; ii++){
                double fit_k1 = (m.prim(IV1, kk, jj, ii) - m.prim(IV1, kk, jj, ii - 1)) / (m.x1v(ii) - m.x1v(ii - 1));
                double fit_b1 = m.prim(IV1, kk, jj, ii) - fit_k1 * m.x1v(ii);
                v1bound1(kk, jj, ii) = fit_k1 * m.x1f(ii) + fit_b1;
            }
        }
    }

    // x2-boundaries
    if (m.dim > 1){
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int jj = m.x2s; jj < m.x2l + 1; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double fit_k2 = (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk, jj - 1, ii)) / (m.x2v(jj) - m.x2v(jj - 1));
                    double fit_b2 = m.prim(IV2, kk, jj, ii) - fit_k2 * m.x2v(jj);
                    v2bound2(kk, jj, ii) = fit_k2 * m.x2f(jj) + fit_b2;
                }
            }
        }
    }
    else{
        v2bound2.set_uniform(0);
    }

    // x3-boundaries
    if (m.dim > 2){
        for (int kk = m.x3s; kk < m.x3l + 1; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double fit_k3 = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk - 1, jj, ii)) / (m.x3v(kk) - m.x3v(kk - 1));
                    double fit_b3 = m.prim(IV3, kk, jj, ii) - fit_k3 * m.x3v(kk);
                    v3bound3(kk, jj, ii) = fit_k3 * m.x3f(kk) + fit_b3;
                }
            }
        }
    }
    else{
        v3bound3.set_uniform(0);
    }

    // Step 2: Use the boundary values to calculate gradient in active domain
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                #if defined(CARTESIAN_COORD)
                gradV(kk, jj, ii) =  (v1bound1(kk, jj, ii) - v1bound1(kk, jj, ii - 1)) / m.dx1p(kk, jj, ii)
                                      + (v2bound2(kk, jj, ii) - v2bound2(kk, jj - 1, ii)) / m.dx2p(kk, jj, ii)
                                      + (v3bound3(kk, jj, ii) - v3bound3(kk - 1, jj, ii)) / m.dx3p(kk, jj, ii);
                #elif defined(SPHERICAL_POLAR_COORD)
                gradV(kk, jj, ii) =  (m.rsq(ii) * v1bound1(kk, jj, ii) - m.rsq(ii - 1) * v1bound1(kk, jj, ii - 1)) / (m.x1v(ii) * m.x1v(ii) * m.dx1(ii))
                                      + (m.geo_sm(jj) * v2bound2(kk, jj, ii) - m.geo_sm(jj - 1) * v2bound2(kk, jj - 1, ii)) / (m.x1v(ii) * std::sin(m.x2v(jj)) * m.dx2(jj))
                                      + (v3bound3(kk, jj, ii) - v3bound3(kk - 1, jj, ii)) / (m.x1v(ii) * std::sin(m.x2v(jj)) * m.dx3(kk));
                #endif // COORD
                /*
                gradV(0, kk, jj, ii) = (v1bound1(kk, jj, ii) - v1bound1(kk, jj, ii - 1)) / m.dx1p(kk, jj, ii)
                                     + (v2bound1(kk, jj, ii) - v2bound1(kk, jj, ii - 1)) / m.dx2p(kk, jj, ii)
                                     + (v3bound1(kk, jj, ii) - v3bound1(kk, jj, ii - 1)) / m.dx3p(kk, jj, ii);
                gradV(1, kk, jj, ii) = (v1bound2(kk, jj, ii) - v1bound2(kk, jj - 1, ii)) / m.dx1p(kk, jj, ii)
                                     + (v2bound2(kk, jj, ii) - v2bound2(kk, jj - 1, ii)) / m.dx2p(kk, jj, ii)
                                     + (v3bound2(kk, jj, ii) - v3bound2(kk, jj - 1, ii)) / m.dx3p(kk, jj, ii);
                gradV(2, kk, jj, ii) = (v1bound3(kk, jj, ii) - v1bound3(kk - 1, jj, ii)) / m.dx1p(kk, jj, ii)
                                     + (v2bound3(kk, jj, ii) - v2bound3(kk - 1, jj, ii)) / m.dx2p(kk, jj, ii)
                                     + (v3bound3(kk, jj, ii) - v3bound3(kk - 1, jj, ii)) / m.dx3p(kk, jj, ii);
                */
            }
        }
    }
    // step 3: Copy the boundary of the active domain into the first ghost zone
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                gradV(kk, jj, m.x1s - 1) = gradV(kk, jj, m.x1s);
            }
        }
    if (m.dim > 1){
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                gradV(kk, m.x2s - 1, ii) = gradV(kk, m.x2s, ii);
            }
        }
    }
    if (m.dim > 2){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                gradV(m.x3s - 1, jj, ii) = gradV(m.x3s, jj, ii);
            }
        }
    }
}


void CalcVisT11(mesh &m, BootesArray<double> &VisT11, BootesArray<double> &gradV){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l + 1; ii++){
                VisT11(kk, jj, ii) = 2 * (m.prim(IV1, kk, jj, ii) - m.prim(IV1, kk, jj, ii - 1)) / m.dx1p(kk, jj, ii) - 2. / 3. * 0.5 * (gradV(kk, jj, ii) + gradV(kk, jj, ii - 1));
            }
        }
    }
}


void CalcVisT12_1only(mesh &m, BootesArray<double> &VisT12){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l + 1; ii++){
                VisT12(kk, jj, ii) = (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk, jj, ii - 1)) / m.dx1p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                // Physically, T12 = flux of M2 in direction of M1, so interpolate along direction 1.
                double fit_k2 = (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk, jj, ii - 1)) / (m.x1v(ii) - m.x1v(ii - 1));
                double fit_b2 = m.prim(IV2, kk, jj, ii) - fit_k2 * m.x1v(ii);
                double v2atbound = fit_k2 * m.x1f(ii) + fit_b2;
                VisT12(kk, jj, ii) -= v2atbound / m.x1f(ii);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT12(mesh &m, BootesArray<double> &VisT12){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l + 1; ii++){
                VisT12(kk, jj, ii) = (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk, jj, ii - 1)) / m.dx1p(kk, jj, ii)
                                   + (m.prim(IV1, kk, jj, ii) - m.prim(IV1, kk, jj - 1, ii)) / m.dx2p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                double fit_k2 = (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk, jj, ii - 1)) / (m.x1v(ii) - m.x1v(ii - 1));
                double fit_b2 = m.prim(IV2, kk, jj, ii) - fit_k2 * m.x1v(ii);
                double v2atbound = fit_k2 * m.x1f(ii) + fit_b2;
                VisT12(kk, jj, ii) -= v2atbound / m.x1f(ii);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT13_1only(mesh &m, BootesArray<double> &VisT13){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l + 1; ii++){
                VisT13(kk, jj, ii) = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj, ii - 1)) / m.dx1p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                double fit_k3 = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj, ii - 1)) / (m.x1v(ii) - m.x1v(ii - 1));
                double fit_b3 = m.prim(IV3, kk, jj, ii) - fit_k3 * m.x1v(ii);
                double v3atbound = fit_k3 * m.x1f(ii) + fit_b3;
                VisT13(kk, jj, ii) -= v3atbound / m.x1f(ii);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT13(mesh &m, BootesArray<double> &VisT13){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l + 1; ii++){
                VisT13(kk, jj, ii) = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj, ii - 1)) / m.dx1p(kk, jj, ii)
                                   + (m.prim(IV1, kk, jj, ii) - m.prim(IV1, kk - 1, jj, ii)) / m.dx3p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                double fit_k3       = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj, ii - 1)) / (m.x1v(ii) - m.x1v(ii - 1));
                double fit_b3       = m.prim(IV3, kk, jj, ii) - fit_k3 * m.x1v(ii);
                double v3atbound    = fit_k3 * m.x1f(ii) + fit_b3;
                VisT13(kk, jj, ii) -= v3atbound / m.x1f(ii);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT21(mesh &m, BootesArray<double> &VisT21){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l + 1; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                VisT21(kk, jj, ii) = (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk, jj, ii - 1)) / m.dx1p(kk, jj, ii)
                                   + (m.prim(IV1, kk, jj, ii) - m.prim(IV1, kk, jj - 1, ii)) / m.dx2p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                // Physically, T21 = flux of M1 in direction 2, so interpolate *the quantity needed* along direction 2
                double fit_k2 = (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk, jj - 1, ii)) / (m.x2v(jj) - m.x2v(jj - 1));
                double fit_b2 = m.prim(IV2, kk, jj, ii) - fit_k2 * m.x2v(jj);
                double v2atbound = fit_k2 * m.x2f(jj) + fit_b2;
                VisT21(kk, jj, ii) -= v2atbound / m.x1v(ii);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT22(mesh &m, BootesArray<double> &VisT22, BootesArray<double> &gradV){
    //VisT22.set_uniform(0);
    //return;
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l + 1; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                VisT22(kk, jj, ii) = 2 * (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk, jj - 1, ii)) / m.dx2p(kk, jj, ii) - 2. / 3. * 0.5 * (gradV(kk, jj - 1, ii) + gradV(kk, jj, ii));
                //#ifdef SPHERICAL_POLAR_COORD
                //double fit_k1 = (m.prim(IV1, kk, jj, ii) - m.prim(IV1, kk, jj - 1, ii)) / (m.x2v(jj) - m.x2v(jj - 1));
                //double fit_b1 = m.prim(IV1, kk, jj, ii) - fit_k1 * m.x2v(jj);
                //double v1atbound = fit_k1 * m.x2f(jj) + fit_b1;

                //VisT22(kk, jj, ii) += 2 * v1atbound / m.x1v(ii);
                //#endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT23_2only(mesh &m, BootesArray<double> &VisT23){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l + 1; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                VisT23(kk, jj, ii) = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj - 1, ii)) / m.dx2p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                double fit_k3 = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj - 1, ii)) / (m.x2v(jj) - m.x2v(jj - 1));
                double fit_b3 = m.prim(IV3, kk, jj, ii) - fit_k3 * m.x2v(jj);
                double v3atbound = fit_k3 * m.x2f(jj) + fit_b3;

                // TODO: the cot term should be cot(x2f(jj)). the one used now is cot(x2v(jj));
                VisT23(kk, jj, ii) -= v3atbound / m.x1v(ii) * m.geo_cot(jj);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT23(mesh &m, BootesArray<double> &VisT23){
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l + 1; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                VisT23(kk, jj, ii) = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj - 1, ii)) / m.dx2p(kk, jj, ii)
                                   + (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk - 1, jj, ii)) / m.dx3p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                double fit_k3 = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj - 1, ii)) / (m.x2v(jj) - m.x2v(jj - 1));
                double fit_b3 = m.prim(IV3, kk, jj, ii) - fit_k3 * m.x2v(jj);
                double v3atbound = fit_k3 * m.x2f(jj) + fit_b3;

                // TODO: the cot term should be cot(x2f(jj)). the one used now is cot(x2v(jj));
                VisT23(kk, jj, ii) -= v3atbound / m.x1v(ii) * m.geo_cot(jj);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT31(mesh &m, BootesArray<double> &VisT31){
    for (int kk = m.x3s; kk < m.x3l + 1; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                VisT31(kk, jj, ii) = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj, ii - 1)) / m.dx1p(kk, jj, ii)
                                   + (m.prim(IV1, kk, jj, ii) - m.prim(IV1, kk - 1, jj, ii)) / m.dx3p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                double fit_k3       = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk - 1, jj, ii)) / (m.x3v(kk) - m.x3v(kk - 1));
                double fit_b3       = m.prim(IV3, kk, jj, ii) - fit_k3 * m.x3v(kk);
                double v3atbound    = fit_k3 * m.x3f(kk) + fit_b3;
                VisT31(kk, jj, ii) -= v3atbound / m.x1v(ii);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT32(mesh &m, BootesArray<double> &VisT32){
    for (int kk = m.x3s; kk < m.x3l + 1; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                VisT32(kk, jj, ii) = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk, jj - 1, ii)) / m.dx2p(kk, jj, ii)
                                   + (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk - 1, jj, ii)) / m.dx3p(kk, jj, ii);
                #ifdef SPHERICAL_POLAR_COORD
                double fit_k3 = (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk - 1, jj, ii)) / (m.x3v(kk) - m.x3v(kk - 1));
                double fit_b3 = m.prim(IV3, kk, jj, ii) - fit_k3 * m.x3v(kk);
                double v3atbound = fit_k3 * m.x3f(kk) + fit_b3;

                VisT32(kk, jj, ii) -= v3atbound / m.x1v(ii) * m.geo_cot(jj);
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void CalcVisT33(mesh &m, BootesArray<double> &VisT33, BootesArray<double> &gradV){
    for (int kk = m.x3s; kk < m.x3l + 1; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                VisT33(kk, jj, ii) = 2 * (m.prim(IV3, kk, jj, ii) - m.prim(IV3, kk - 1, jj, ii)) / m.dx3p(kk, jj, ii) - 2. / 3. * 0.5 * (gradV(kk - 1, jj, ii) + gradV(kk, jj, ii));
                #ifdef SPHERICAL_POLAR_COORD
                double fit_k1 = (m.prim(IV1, kk, jj, ii) - m.prim(IV1, kk - 1, jj, ii)) / (m.x3v(kk) - m.x3v(kk - 1));
                double fit_b1 = m.prim(IV1, kk, jj, ii) - fit_k1 * m.x3v(kk);
                double v1atbound = fit_k1 * m.x3f(kk) + fit_b1;
                double fit_k2 = (m.prim(IV2, kk, jj, ii) - m.prim(IV2, kk - 1, jj, ii)) / (m.x3v(kk) - m.x3v(kk - 1));
                double fit_b2 = m.prim(IV2, kk, jj, ii) - fit_k2 * m.x3v(kk);
                double v2atbound = fit_k2 * m.x3f(kk) + fit_b2;

                VisT33(kk, jj, ii) -= (2. * v2atbound / m.x1v(ii) + v2atbound / m.x1v(ii) * m.geo_cot(jj));
                #endif // SPHERICAL_POLAR_COORD
            }
        }
    }
}


void apply_viscous_flux(mesh &m, double &dt, BootesArray<double> &fcons, BootesArray<double> &nu_vis){
    // step 1: make arrays storing grad*V and the viscous stress tensor
    BootesArray<double> gradV;
    BootesArray<double> VisTensor11;
    BootesArray<double> VisTensor12;
    BootesArray<double> VisTensor13;
    BootesArray<double> VisTensor21;
    BootesArray<double> VisTensor22;
    BootesArray<double> VisTensor23;
    BootesArray<double> VisTensor31;
    BootesArray<double> VisTensor32;
    BootesArray<double> VisTensor33;
    gradV.NewBootesArray(m.x3v.shape()[0], m.x2v.shape()[0], m.x1v.shape()[0]);
    // The shapes are different because computation is un-necessary outside the simulation domain.
    VisTensor11.NewBootesArray(m.x3v.shape()[0], m.x2v.shape()[0], m.x1f.shape()[0]);       // in x1 direction, flux of M1
    VisTensor12.NewBootesArray(m.x3v.shape()[0], m.x2v.shape()[0], m.x1f.shape()[0]);       // in x1 direction, flux of M2
    VisTensor13.NewBootesArray(m.x3v.shape()[0], m.x2v.shape()[0], m.x1f.shape()[0]);       // in x1 direction, flux of M3
    VisTensor21.NewBootesArray(m.x3v.shape()[0], m.x2f.shape()[0], m.x1v.shape()[0]);       // in x2 direction, flux of M1
    VisTensor22.NewBootesArray(m.x3v.shape()[0], m.x2f.shape()[0], m.x1v.shape()[0]);       // in x2 direction, flux of M2
    VisTensor23.NewBootesArray(m.x3v.shape()[0], m.x2f.shape()[0], m.x1v.shape()[0]);       // in x2 direction, flux of M3
    VisTensor31.NewBootesArray(m.x3f.shape()[0], m.x2v.shape()[0], m.x1v.shape()[0]);       // in x3 direction, flux of M1
    VisTensor32.NewBootesArray(m.x3f.shape()[0], m.x2v.shape()[0], m.x1v.shape()[0]);       // in x3 direction, flux of M2
    VisTensor33.NewBootesArray(m.x3f.shape()[0], m.x2v.shape()[0], m.x1v.shape()[0]);       // in x3 direction, flux of M3
    // step 2: fill in values for grad*V
    calcgradV(m, gradV);
    // step 3: fill in values for the stress tensor
    if (m.dim == 1){
        CalcVisT11(m, VisTensor11, gradV);
        CalcVisT12_1only(m, VisTensor12);
        CalcVisT13_1only(m, VisTensor13);
    }
    if (m.dim == 2){
        CalcVisT11(m, VisTensor11, gradV);
        CalcVisT12(m, VisTensor12);
        CalcVisT13_1only(m, VisTensor13);
        CalcVisT21(m, VisTensor21);
        CalcVisT22(m, VisTensor22, gradV);
        CalcVisT23_2only(m, VisTensor23);
    }
    if (m.dim == 3){
        CalcVisT11(m, VisTensor11, gradV);
        CalcVisT12(m, VisTensor12);
        CalcVisT13(m, VisTensor13);
        CalcVisT21(m, VisTensor21);
        CalcVisT22(m, VisTensor22, gradV);
        CalcVisT23(m, VisTensor23);
        CalcVisT31(m, VisTensor31);
        CalcVisT32(m, VisTensor32);
        CalcVisT33(m, VisTensor33, gradV);
    }
    std::cout << "finish calc VIS" << std::endl << std::flush;
    // step 4: use the stress tensor values to update flux of the conservatives
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l + 1; ii++){
                int kkf = kk - m.x3s;
                int jjf = jj - m.x2s;
                int iif = ii - m.x1s;
                fcons(IM1, 0, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor11(kk, jj, ii);
                fcons(IM2, 0, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor12(kk, jj, ii);
                fcons(IM3, 0, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor13(kk, jj, ii);
            }
        }
    }
    std::cout << "finish 1" << std::endl << std::flush;
    if (m.dim > 1){
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int jj = m.x2s; jj < m.x2l + 1; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    int kkf = kk - m.x3s;
                    int jjf = jj - m.x2s;
                    int iif = ii - m.x1s;
                    fcons(IM1, 1, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor21(kk, jj, ii);
                    fcons(IM2, 1, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor22(kk, jj, ii);
                    fcons(IM3, 1, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor23(kk, jj, ii);
                }
            }
        }
    }
    std::cout << "finish 2" << std::endl << std::flush;
    if (m.dim > 2){
        for (int kk = m.x3s; kk < m.x3l + 1; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    int kkf = kk - m.x3s;
                    int jjf = jj - m.x2s;
                    int iif = ii - m.x1s;
                    fcons(IM1, 2, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor31(kk, jj, ii);
                    fcons(IM2, 2, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor32(kk, jj, ii);
                    fcons(IM3, 2, kkf, jjf, iif) += m.prim(IDN, kk, jj, ii) * nu_vis(kk, jj, ii) * VisTensor33(kk, jj, ii);
                }
            }
        }
    }
    std::cout << "finish 3" << std::endl << std::flush;
}


