#ifndef MESH_HPP_
#define MESH_HPP_
#include <string>
#include <memory>
#include <iostream>
#include "BootesArray.hpp"
#include "util.hpp"
#include "eos/eos.hpp"
#include "boundary_condition/standard_bc.hpp"

int NUMCONS = 5;
int NUMPRIM = 5;
enum ConsIndex:int{IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum PrimIndex:int{IDP=0, IV1=1, IV2=2, IV3=3, IPN=4};

class mesh{
    public:
        /** grid **/
        int dim;
        BootesArray<float> x1v;       // cell center (1D array)
        BootesArray<float> x2v;
        BootesArray<float> x3v;
        BootesArray<float> x1f;       // cell face (1D array)
        BootesArray<float> x2f;
        BootesArray<float> x3f;
        BootesArray<float> dx1;       // delta x (1D array), coordinate unit
        BootesArray<float> dx2;
        BootesArray<float> dx3;
        BootesArray<float> dx1p;      // delta x (1D array), physical unit
        BootesArray<float> dx2p;
        BootesArray<float> dx3p;

        // for spherical polar coordinate, since dxi is not sufficient.
        BootesArray<float> vol;          // volume of each cell
        BootesArray<float> f1a;          // face size in x1 (3D array), size (Nx3, Nx2, Nx1 + 1), Nxi = nxi + 2 * ngi
        BootesArray<float> f2a;          //
        BootesArray<float> f3a;          //
        BootesArray<float> one_orgeo;    // = 0.5 * (rp^2 - rm^2) / (1/3 * rp^3 - rm^3) used in geometry terms
        BootesArray<float> rV;           // = (rm + rp) * 1/3 * (rp^3 - rm^3)
        BootesArray<float> geo_cot;         // = (sin(tp) - sin(tm)) / abs(cos(tp) - cos(tm))
        BootesArray<float> geo_sm;         // = (sin(tp) - sin(tm)) / abs(cos(tp) - cos(tm))
        BootesArray<float> geo_sp;         // = (sin(tp) - sin(tm)) / abs(cos(tp) - cos(tm))

        int x1s, x2s, x3s;                     // start index of active domain
        int x1l, x2l, x3l;                     // end index of active domain
        int nx1, nx2, nx3;                     // number of active zones in each direction
        int ng1, ng2, ng3;                     // number of ghost zones in each direction, implement for 2D and 1D simulation

        float hydro_gamma;

        /** cons **/
        BootesArray<float> cons;           // 4D (5, z, y, x)

        /** prim **/
        BootesArray<float> prim;            // 4D (4, z, y, x)

        /** functions **/
        void prim_to_cons();
        void cons_to_prim();

        /** setup grid functions **/
        void SetupCartesian(int dimension,
                            float x1min, float x1max, int numx1, int ngh1,
                            float x2min, float x2max, int numx2, int ngh2,
                            float x3min, float x3max, int numx3, int ngh3);
        void SetupSphericalPolar(int dimension,
                                 float x1min, float x1max, int numx1, float ratio1, int ngh1,
                                 float x2min, float x2max, int numx2,               int ngh2,
                                 float x3min, float x3max, int numx3,               int ngh3);
};



void mesh::SetupCartesian(int dimension,
                          float x1min, float x1max, int numx1, int ngh1,
                          float x2min, float x2max, int numx2, int ngh2,
                          float x3min, float x3max, int numx3, int ngh3){
    dim = dimension;
    ng1 = ngh1;  ng2 = ngh2;  ng3 = ngh3;
    nx1 = numx1; nx2 = numx2; nx3 = numx3;
    x1s = ng1; x1l = numx1 + ng1;
    x2s = ng2; x2l = numx2 + ng2;
    x3s = ng3; x3l = numx3 + ng3;
    float dx1_num = (x1max - x1min) / numx1;
    float dx2_num = (x2max - x2min) / numx2;
    float dx3_num = (x3max - x3min) / numx3;

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
    dx1p = dx1;
    dx2p = dx2;
    dx3p = dx3;

    // allocate memory for hydro quantities
    cons.NewBootesArray(NUMCONS, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    prim.NewBootesArray(NUMPRIM, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
}


void mesh::SetupSphericalPolar(int dimension,
                               float x1min, float x1max, int numx1, float ratio1, int ngh1,
                               float x2min, float x2max, int numx2,               int ngh2,
                               float x3min, float x3max, int numx3,               int ngh3){
    // 1 - r; 2 - theta; 3 - phi
    dim = dimension;
    ng1 = ngh1;  ng2 = ngh2;  ng3 = ngh3;
    nx1 = numx1; nx2 = numx2; nx3 = numx3;
    x1s = ng1; x1l = numx1 + ng1;
    x2s = ng2; x2l = numx2 + ng2;
    x3s = ng3; x3l = numx3 + ng3;
    float dx2_num = (x2max - x2min) / numx2;
    float dx3_num = (x3max - x3min) / numx3;
    float mx1 = (x1max - x1min) * (1 - ratio1) / (1 - pow(ratio1, numx1));

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

    // step 3: setup cell size values
    dx1.NewBootesArray(x1f.shape()[0] - 1);
    dx2.NewBootesArray(x2f.shape()[0] - 1);
    dx3.NewBootesArray(x3f.shape()[0] - 1);
    for (int ii = 0; ii < dx1.shape()[0]; ii ++) { dx1(ii) = x1f(ii + 1) - x1f(ii); }
    for (int ii = 0; ii < dx2.shape()[0]; ii ++) { dx2(ii) = x2f(ii + 1) - x2f(ii); }
    for (int ii = 0; ii < dx3.shape()[0]; ii ++) { dx3(ii) = x3f(ii + 1) - x3f(ii); }
    dx1p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    dx2p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    dx3p.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    // TODO: replace x1v(ii) with higher order ones.
    for (int kk = 0; kk < x3v.shape()[0]; kk ++){
        for (int jj = 0; jj < x2v.shape()[0]; jj ++){
            for (int ii = 0; ii < x1v.shape()[0]; ii ++){
                dx1p(kk, jj, ii) = dx1(ii);
                dx2p(kk, jj, ii) = dx2(jj) * x1v(ii);
                dx3p(kk, jj, ii) = dx3(kk) * 1.0 / x1v(ii) * -1 * (cos(x2f(jj)) - cos(x2f(jj + 1))) / (dx2(jj));
            }
        }
    }

    // step 3.1: setup cell size value in physical units

    // step 4: setup face sizes
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

    // step 5: setup volume size of each cell
    vol.NewBootesArray(x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    for (int kk = 0; kk < x3v.shape()[0]; kk ++){
        for (int jj = 0; jj < x2v.shape()[0]; jj ++){
            for (int ii = 0; ii < x1v.shape()[0]; ii ++){
                vol(kk, jj, ii) = 1. / 3. * (pow(x1f(ii + 1), 3) - pow(x1f(ii), 3)) * (cos(x2f(jj)) - cos(x2f(jj + 1))) * (x3f(kk + 1) - x3f(kk));
            }
        }
    }

    // step 6: setup other geometric terms
    // step 6.1: radial direction
    one_orgeo.NewBootesArray(x1v.shape()[0]);
    rV.NewBootesArray(x1v.shape()[0]);
    for (int ii = 0; ii < x1v.shape()[0]; ii ++ ){
        float rp = x1f(ii + 1);
        float rm = x1f(ii);
        one_orgeo(ii) = 0.5 * (rp * rp - rm * rm) / (1./3. * (rp * rp * rp - rm * rm * rm));
        rV(ii) = (rm + rp) * 1. / 3. * (rp * rp * rp - rm * rm * rm);
    }
    // step 6.2: theta direction
    geo_cot.NewBootesArray(x2v.shape()[0]);
    geo_sm.NewBootesArray(x2v.shape()[0]);
    geo_sp.NewBootesArray(x2v.shape()[0]);
    for (int jj = 0; jj < geo_cot.shape()[0]; jj ++){
        float sm = std::sin(x2f(jj));
        float sp = std::sin(x2f(jj+1));
        float cm = std::cos(x2f(jj));
        float cp = std::cos(x2f(jj+1));
        geo_sm(jj) = sm;
        geo_sp(jj) = sp;
        geo_cot(jj) = (sp - sm) / abs(cm - cp);
    }


    // step 7: allocate memory for hydro quantities
    cons.NewBootesArray(NUMCONS, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
    prim.NewBootesArray(NUMPRIM, x3v.shape()[0], x2v.shape()[0], x1v.shape()[0]);
}


void mesh::prim_to_cons(){
    ;
}


void mesh::cons_to_prim(){
    #pragma omp parallel for collapse (3) schedule (static)
    for (int kk = x3s; kk < x3l ; kk++){
        for (int jj = x2s; jj < x2l; jj++){
            for (int ii = x1s; ii < x1l; ii++){
                float v1 = cons(IM1, kk, jj, ii) / cons(IDN, kk, jj, ii);
                float v2 = cons(IM2, kk, jj, ii) / cons(IDN, kk, jj, ii);
                float v3 = cons(IM3, kk, jj, ii) / cons(IDN, kk, jj, ii);
                prim(IDN, kk, jj, ii) = cons(IDN, kk, jj, ii);
                prim(IV1, kk, jj, ii) = v1;
                prim(IV2, kk, jj, ii) = v2;
                prim(IV3, kk, jj, ii) = v3;
                prim(IPN, kk, jj, ii) = pres(cons(IDN, kk, jj, ii), cons(IEN, kk, jj, ii),
                                             cons(IM1, kk, jj, ii), cons(IM2, kk, jj, ii), cons(IM3, kk, jj, ii), hydro_gamma);
            }
        }
    }
}


#endif
