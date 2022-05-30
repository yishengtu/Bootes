#include "mesh.hpp"
#include "../index_def.hpp"
#include "../gravity/gravity.hpp"
#include "../../defs.hpp"

void mesh::SetupCartesian(int dimension,
                          double x1min, double x1max, int numx1, int ngh1,
                          double x2min, double x2max, int numx2, int ngh2,
                          double x3min, double x3max, int numx3, int ngh3){
    dim = dimension;
    ng1 = ngh1;  ng2 = ngh2;  ng3 = ngh3;
    nx1 = numx1; nx2 = numx2; nx3 = numx3;
    x1s = ng1; x1l = numx1 + ng1;
    x2s = ng2; x2l = numx2 + ng2;
    x3s = ng3; x3l = numx3 + ng3;
    double dx1_num = (x1max - x1min) / numx1;
    double dx2_num = (x2max - x2min) / numx2;
    double dx3_num = (x3max - x3min) / numx3;

    x1f = linspace(x1min - ng1 * dx1_num, x1max + ng1 * dx1_num, numx1 + 2 * ng1 + 1, true);
    x2f = linspace(x2min - ng2 * dx2_num, x2max + ng2 * dx2_num, numx2 + 2 * ng2 + 1, true);
    x3f = linspace(x3min - ng3 * dx3_num, x3max + ng3 * dx3_num, numx3 + 2 * ng3 + 1, true);

    x1v.NewBootesArray(x1f.shape()[0] - 1);
    x2v.NewBootesArray(x2f.shape()[0] - 1);
    x3v.NewBootesArray(x3f.shape()[0] - 1);

    for (int ii = 0; ii < x1v.shape()[0]; ii ++) { x1v(ii) = (x1f(ii) + x1f(ii + 1)) / 2.0; }
    for (int ii = 0; ii < x2v.shape()[0]; ii ++) { x2v(ii) = (x2f(ii) + x2f(ii + 1)) / 2.0; }
    for (int ii = 0; ii < x3v.shape()[0]; ii ++) { x3v(ii) = (x3f(ii) + x3f(ii + 1)) / 2.0; }
    dx1.NewBootesArray(x1f.shape()[0] - 1);
    dx2.NewBootesArray(x2f.shape()[0] - 1);
    dx3.NewBootesArray(x3f.shape()[0] - 1);
    for (int ii = 0; ii < dx1.shape()[0]; ii ++) { dx1(ii) = x1f(ii + 1) - x1f(ii); }
    for (int ii = 0; ii < dx2.shape()[0]; ii ++) { dx2(ii) = x2f(ii + 1) - x2f(ii); }
    for (int ii = 0; ii < dx3.shape()[0]; ii ++) { dx3(ii) = x3f(ii + 1) - x3f(ii); }
    dx1p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    dx2p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    dx3p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    // step 5.2: setup cell size value in physical units
    // TODO: use higher order terms
    for (int kk = 0; kk < x3v.shape()[0]; kk ++){
        for (int jj = 0; jj < x2v.shape()[0]; jj ++){
            for (int ii = 0; ii < x1v.shape()[0]; ii ++){
                dx1p(kk, jj, ii) = dx1(ii);
                dx2p(kk, jj, ii) = dx2(jj);
                dx3p(kk, jj, ii) = dx3(kk);
            }
        }
    }

    // allocate memory for hydro quantities
    cons.NewBootesArray(NUMCONS, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    prim.NewBootesArray(NUMPRIM, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);

    // allocate memory for other necessary fields
    #if defined (ENABLE_GRAVITY)
        grav->setup_Phimesh(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    #endif // defined
}


void mesh::SetupSphericalPolar(int dimension,
                               double x1min, double x1max, int numx1, double ratio1, int ngh1,
                               double x2min, double x2max, int numx2,                int ngh2,
                               double x3min, double x3max, int numx3,                int ngh3){
    // 1 - r; 2 - theta; 3 - phi
    dim = dimension;
    ng1 = ngh1;  ng2 = ngh2;  ng3 = ngh3;
    nx1 = numx1; nx2 = numx2; nx3 = numx3;
    x1s = ng1; x1l = numx1 + ng1;
    x2s = ng2; x2l = numx2 + ng2;
    x3s = ng3; x3l = numx3 + ng3;
    double dx2_num = (x2max - x2min) / numx2;
    double dx3_num = (x3max - x3min) / numx3;
    double mx1 = (x1max - x1min) * (1 - ratio1) / (1 - pow(ratio1, numx1));

    // step 1: setup boundary locations
    x1f.NewBootesArray(nx1 + 2 * ng1 + 1);
    x1f(0) = x1min;
    for (int ii = 1; ii < x1f.shape()[0]; ii++){
        x1f(ii) = x1f(ii - 1) + mx1 * pow(ratio1, ii);
    }
    x2f = linspace(x2min - ng2 * dx2_num, x2max + ng2 * dx2_num, numx2 + 2 * ng2 + 1, true);
    x3f = linspace(x3min - ng3 * dx3_num, x3max + ng3 * dx3_num, numx3 + 2 * ng3 + 1, true);

    // step 2: setup cell center locations
    x1v.NewBootesArray(x1f.shape()[0] - 1);
    x2v.NewBootesArray(x2f.shape()[0] - 1);
    x3v.NewBootesArray(x3f.shape()[0] - 1);

    for (int ii = 0; ii < x1v.shape()[0]; ii ++) {
        x1v(ii) = 3.0 / 4.0 * (pow(x1f(ii + 1), 4) - pow(x1f(ii), 4)) / (pow(x1f(ii + 1), 3) - pow(x1f(ii), 3));
    }
    for (int jj = 0; jj < x2v.shape()[0]; jj ++) {
        x2v(jj) = (x2f(jj) * cos(x2f(jj)) - x2f(jj + 1) * cos(x2f(jj + 1)) - sin(x2f(jj)) + sin(x2f(jj + 1))) / (cos(x2f(jj)) - cos(x2f(jj + 1)));
    }
    for (int ii = 0; ii < x3v.shape()[0]; ii ++) {
        x3v(ii) = (x3f(ii) + x3f(ii + 1)) / 2.0;
    }


    // step 3: setup face sizes
    f1a.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1f.shape()[0]);
    f2a.NewBootesArray(x3v.shape()[0], x2f.shape()[0], x1v.shape()[0]);
    f3a.NewBootesArray(x3f.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    #pragma omp parallel for collapse (3)
    for (int kk = 0; kk < x3f.shape()[0]; kk ++){
        for (int jj = 0; jj < x2f.shape()[0]; jj ++){
            for (int ii = 0; ii < x1f.shape()[0]; ii ++){
                if (jj != x2f.shape()[0] - 1 && kk != x3f.shape()[0] - 1){
                    f1a(kk, jj, ii) = pow(x1f(ii), 2) * (cos(x2f(jj)) - cos(x2f(jj + 1))) * (x3f(kk + 1) - x3f(kk));
                }
                if (ii != x1f.shape()[0] - 1 && kk != x3f.shape()[0] - 1){
                    f2a(kk, jj, ii) = 0.5 * (pow(x1f(ii + 1), 2) - pow(x1f(ii), 2)) * sin(x2f(jj)) * (x3f(kk + 1) - x3f(kk));
                }
                if (ii != x1f.shape()[0] - 1 && jj != x2f.shape()[0] - 1){
                    f3a(kk, jj, ii) = 0.5 * (pow(x1f(ii + 1), 2) - pow(x1f(ii), 2)) * (x2f(jj + 1) - x2f(jj));
                }
            }
        }
    }

    // step 4: setup volume size of each cell
    vol.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    for (int kk = 0; kk < x3v.shape()[0]; kk ++){
        for (int jj = 0; jj < x2v.shape()[0]; jj ++){
            for (int ii = 0; ii < x1v.shape()[0]; ii ++){
                vol(kk, jj, ii) = 1. / 3. * (pow(x1f(ii + 1), 3) - pow(x1f(ii), 3)) * (cos(x2f(jj)) - cos(x2f(jj + 1))) * (x3f(kk + 1) - x3f(kk));
            }
        }
    }

    // step 5: setup cell size values
    // step 5.1: setup cell size value in coordinate units
    dx1.NewBootesArray(x1f.shape()[0] - 1);
    dx2.NewBootesArray(x2f.shape()[0] - 1);
    dx3.NewBootesArray(x3f.shape()[0] - 1);
    for (int ii = 0; ii < dx1.shape()[0]; ii ++) { dx1(ii) = x1f(ii + 1) - x1f(ii); }
    for (int ii = 0; ii < dx2.shape()[0]; ii ++) { dx2(ii) = x2f(ii + 1) - x2f(ii); }
    for (int ii = 0; ii < dx3.shape()[0]; ii ++) { dx3(ii) = x3f(ii + 1) - x3f(ii); }
    dx1p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    dx2p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    dx3p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    // step 5.2: setup cell size value in physical units
    // TODO: use higher order terms
    for (int kk = 0; kk < x3v.shape()[0]; kk ++){
        for (int jj = 0; jj < x2v.shape()[0]; jj ++){
            for (int ii = 0; ii < x1v.shape()[0]; ii ++){
                dx1p(kk, jj, ii) = dx1(ii);
                dx2p(kk, jj, ii) = dx2(jj) * x1v(ii);
                dx3p(kk, jj, ii) = dx3(kk) * x1v(ii) * (cos(x2f(jj)) - cos(x2f(jj + 1))) / (dx2(jj));
            }
        }
    }


    // step 6: setup other geometric terms
    // step 6.1: radial direction
    one_orgeo.NewBootesArray(x1v.shape()[0]);
    rV.NewBootesArray(x1v.shape()[0]);
    for (int ii = 0; ii < x1v.shape()[0]; ii ++ ){
        double rp = x1f(ii + 1);
        double rm = x1f(ii);
        one_orgeo(ii) = 0.5 * (rp * rp - rm * rm) / (1./3. * (rp * rp * rp - rm * rm * rm));
        rV(ii) = (rm + rp) * 1. / 3. * (rp * rp * rp - rm * rm * rm);
    }
    rsq.NewBootesArray(x1f.shape()[0]);
    for (int ii = 0; ii < x1f.shape()[0]; ii ++){
        rsq(ii) = x1f(ii) * x1f(ii);
    }
    // step 6.2: theta direction
    geo_cot.NewBootesArray(x2v.shape()[0]);
    geo_sm.NewBootesArray(x2v.shape()[0]);
    geo_sp.NewBootesArray(x2v.shape()[0]);
    for (int jj = 0; jj < geo_cot.shape()[0]; jj ++){
        double sm = std::sin(x2f(jj));
        double sp = std::sin(x2f(jj+1));
        double cm = std::cos(x2f(jj));
        double cp = std::cos(x2f(jj+1));
        geo_sm(jj) = sm;
        geo_sp(jj) = sp;
        geo_cot(jj) = (sp - sm) / (cm - cp);
    }


    // step 7: allocate memory for hydro quantities
    cons.NewBootesArray(NUMCONS, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    prim.NewBootesArray(NUMPRIM, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    // allocate memory for other necessary fields
    #if defined (ENABLE_GRAVITY)
        grav->setup_Phimesh(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    #endif // defined
}


#if defined (ENABLE_DUSTFLUID)
void mesh::setupDustFluidMesh(int NS){
    NUMSPECIES = NS;
    dcons.NewBootesArray(NS, NUMCONS, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    dprim.NewBootesArray(NS, NUMPRIM, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
}
#endif
