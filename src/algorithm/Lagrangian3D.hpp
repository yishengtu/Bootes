#ifndef LAGRANGIAN_3D_HPP_
#define LAGRANGIAN_3D_HPP_
#include <cmath>
#include <string>
#include <hdf5.h>
#include <iostream>
#include <iomanip>
#include "physical_constants.hpp"

#include "BootesArray.hpp"
#include "athena_read.hpp"
#include "util.hpp"


using namespace std;


class CartesianCoord{
    public:
        int nx;
        int ny;
        int nz;
        float xmin, xmax;
        float ymin, ymax;
        float zmin, zmax;

        BootesArray<float> xlist;
        BootesArray<float> ylist;
        BootesArray<float> zlist;

        BootesArray<float> dxlist;
        BootesArray<float> dylist;
        BootesArray<float> dzlist;

        BootesArray<float> rho;
        BootesArray<float> temp;
        BootesArray<float> vx;
        BootesArray<float> vy;
        BootesArray<float> vz;
        BootesArray<float> gravx;
        BootesArray<float> gravy;
        BootesArray<float> gravz;

        BootesArray<float> gradient(BootesArray<float> &y, BootesArray<float> &x){
            BootesArray<float> dydx;
            dydx.NewBootesArray(y.shape()[0]);
            for (int ii = 0; ii < y.shape()[0]; ii++){
                if (ii = 0){
                    dydx(ii) = (y(1) - y(0)) / (x(1) - x(0));
                }
                else if (ii == y.shape()[0] - 1){
                    int last_ind = y.shape()[0] - 1;
                    int sec_last_ind = y.shape()[0] - 2;
                    dydx(ii) = (y(sec_last_ind) - y(last_ind)) / (x(sec_last_ind) - x(last_ind));
                }
                else{
                    dydx(ii) = (y(ii + 1) - y(ii - 1)) / (x(ii + 1) - x(ii - 1));
                }
            }
            return dydx;
        }

        void calc_grav_acceleration(BootesArray<float> &Phi,   BootesArray<float> &xvlist, BootesArray<float> &yvlist, BootesArray<float> &zvlist,
                                    BootesArray<float> &gravx, BootesArray<float> &gravy, BootesArray<float> &gravz){
            /** cx, cy and cz are cell-center values **/
            gravx.NewBootesArray(Phi.shape()[0], Phi.shape()[1], Phi.shape()[2]);
            gravy.NewBootesArray(Phi.shape()[0], Phi.shape()[1], Phi.shape()[2]);
            gravz.NewBootesArray(Phi.shape()[0], Phi.shape()[1], Phi.shape()[2]);

            #pragma omp parallel for collapse (3) schedule(static)
            for (int kk = 0; kk < Phi.shape()[0]; kk++){
                for (int jj = 0; jj < Phi.shape()[1]; jj++){
                    for (int ii = 0; ii < Phi.shape()[2]; ii++){
                        if (ii == 0)                      { gravx(kk, jj, ii) = (Phi(kk, jj, ii + 1) - Phi(kk, jj, ii))     / (xvlist(ii + 1) - xvlist(ii)); }
                        else if (ii == Phi.shape()[2] - 1){ gravx(kk, jj, ii) = (Phi(kk, jj, ii)     - Phi(kk, jj, ii - 1)) / (xvlist(ii)     - xvlist(ii - 1)); }
                        else                              { gravx(kk, jj, ii) = (Phi(kk, jj, ii + 1) - Phi(kk, jj, ii - 1)) / (xvlist(ii + 1) - xvlist(ii - 1)); }

                        if (jj == 0)                      { gravy(kk, jj, ii) = (Phi(kk, jj + 1, ii) - Phi(kk, jj, ii))     / (yvlist(jj + 1) - yvlist(jj)); }
                        else if (jj == Phi.shape()[1] - 1){ gravy(kk, jj, ii) = (Phi(kk, jj, ii)     - Phi(kk, jj - 1, ii)) / (yvlist(jj)     - yvlist(jj - 1)); }
                        else                              { gravy(kk, jj, ii) = (Phi(kk, jj + 1, ii) - Phi(kk, jj - 1, ii)) / (yvlist(jj + 1) - yvlist(jj - 1)); }

                        if (kk == 0)                      { gravz(kk, jj, ii) = (Phi(kk + 1, jj, ii) - Phi(kk, jj, ii))     / (zvlist(kk + 1) - zvlist(kk)); }
                        else if (kk == Phi.shape()[0] - 1){ gravz(kk, jj, ii) = (Phi(kk, jj, ii)     - Phi(kk - 1, jj, ii)) / (zvlist(kk)     - zvlist(kk - 1)); }
                        else                              { gravz(kk, jj, ii) = (Phi(kk + 1, jj, ii) - Phi(kk - 1, jj, ii)) / (zvlist(kk + 1) - zvlist(kk - 1)); }
                    }
                }
            }
            #pragma omp barrier
        }


        void set_hydro_with_time(OutputAthdf &f1, OutputAthdf &f2){
            float time1 = f1.Time;
            float time2 = f2.Time;
            if (time1 == time2){
                cout << "Warning: consecutive frame have identical time!" << endl << flush;
                time2 += 100;
            }
            BootesArray<float> f1_gravx;
            BootesArray<float> f1_gravy;
            BootesArray<float> f1_gravz;
            BootesArray<float> f2_gravx;
            BootesArray<float> f2_gravy;
            BootesArray<float> f2_gravz;

            calc_grav_acceleration(f1.Phi, f1.x1v, f1.x2v, f1.x3v, f1_gravx, f1_gravy, f1_gravz);
            calc_grav_acceleration(f2.Phi, f2.x1v, f2.x2v, f2.x3v, f2_gravx, f2_gravy, f2_gravz);

            set_hydro_two_frames(time1, time2, f1.rho,  f2.rho,  rho);
            set_hydro_two_frames(time1, time2, f1.temp, f2.temp, temp);
            set_hydro_two_frames(time1, time2, f1.vel1, f2.vel1, vx);
            set_hydro_two_frames(time1, time2, f1.vel2, f2.vel2, vy);
            set_hydro_two_frames(time1, time2, f1.vel3, f2.vel3, vz);
            set_hydro_two_frames(time1, time2, f1_gravx, f2_gravx, gravx);
            set_hydro_two_frames(time1, time2, f1_gravy, f2_gravy, gravy);
            set_hydro_two_frames(time1, time2, f1_gravz, f2_gravz, gravz);
        }

        void set_hydro_two_frames(float hydrotime1, float hydrotime2, BootesArray<float> &hydro1, BootesArray<float> &hydro2, BootesArray<float> &mesh_out){
            mesh_out.NewBootesArray(hydro1.shape()[0], hydro1.shape()[1], hydro1.shape()[2], 2);
            float zero = 0.0;
            float hydro_time_diff = hydrotime2 - hydrotime1;
            #pragma omp parallel for collapse (3) schedule(static)
            for (int kk = 0; kk < hydro1.shape()[0]; kk++){
                for (int jj = 0; jj < hydro1.shape()[1]; jj++){
                    for (int ii = 0; ii < hydro1.shape()[2]; ii++){
                        interp_slope(zero, hydro_time_diff, hydro1(kk, jj, ii), hydro2(kk, jj, ii), mesh_out(kk, jj, ii, 0), mesh_out(kk, jj, ii, 1));
                        // cout << hydrotime1 << '\t' << hydrotime2 << '\t' << mesh_out(kk, jj, ii, 0) << '\t' << mesh_out(kk, jj, ii, 1) << endl;
                    }
                }
            }
            #pragma omp barrier
        }

        void interp_slope(float &x1, float &x2, float &y1, float &y2, float &a, float &b){
            a = (y2 - y1) / (x2 - x1);
            b = y1 - a * x1;
        }
};


void TSC(float &x, float &y, float &z,
         BootesArray<float> &xlist, BootesArray<float> &dxlist,
         BootesArray<float> &ylist, BootesArray<float> &dylist,
         BootesArray<float> &zlist, BootesArray<float> &dzlist,
         float w[3][3][3], int& indi, int& indj, int& indk){
    indi = searchsorted(x, xlist) - 1;
    indj = searchsorted(y, ylist) - 1;
    indk = searchsorted(z, zlist) - 1;
    if (indi == -1) { indi = 0; }
    if (indj == -1) { indj = 0; }
    if (indk == -1) { indk = 0; }

    float x1 = (x - xlist(indi)) / dxlist(indi);
    float x2 = (y - ylist(indj)) / dylist(indj);
    float x3 = (z - zlist(indk)) / dzlist(indk);

    float wx1l = 0.5 * (1. - x1) * (1. - x1);
    float wx1r = 0.5 * x1 * x1;
    float wx1m = 1. - (wx1l + wx1r);
    float wx2l = 0.5 * (1. - x2) * (1. - x2);
    float wx2r = 0.5 * x2 * x2;
    float wx2m = 1. - (wx2l + wx2r);
    float wx3l = 0.5 * (1. - x3) * (1. - x3);
    float wx3r = 0.5 * x3 * x3;
    float wx3m = 1. - (wx3l + wx3r);

    float wx1[3] = {wx1l, wx1m, wx1r};
    float wx2[3] = {wx2l, wx2m, wx2r};
    float wx3[3] = {wx3l, wx3m, wx3r};

    for (int kk = 0; kk < 3; kk ++){
        for (int jj = 0; jj < 3; jj ++){
            for (int ii = 0; ii < 3; ii ++){
                w[kk][jj][ii] = wx3[kk] * wx2[jj] * wx1[ii];
            }
        }
    }
}


float cs(float &T){
    return sqrt(kb * T / (hydro_mu * mH));
}


float vth(float &T){
    return sqrt(8. * kb * T / (M_PI * hydro_mu * mH));
}


float apply_TSC_slope_periodicBC(float w[3][3][3], int &i, int &j, int &k, BootesArray<float> &mesh, int &nx, int &ny, int &nz, float &t){
    float summed_quan = 0;
    for (int zi = -1; zi < 2; zi ++){
        int indz = k + zi;
        if (indz > nz - 1){ indz -= nz; }
        else if (indz < 0){ indz += nz; }
        else {;}
        for (int yi = -1; yi < 2; yi ++){
            int indy = j + yi;
            if (indy > ny - 1){ indy -= ny; }
            else if (indy < 0){ indy += ny; }
            else {;}
            for (int xi = -1; xi < 2; xi ++){
                int indx = i + xi;
                if (indx > nx - 1){ indx -= nx; }
                else if (indx < 0){ indx += nx; }
                else {;}

                summed_quan += (mesh(indz, indy, indx, 0) * t + mesh(indz, indy, indx, 1)) * w[xi+1][yi+1][zi+1];
            }
        }
    }
    return summed_quan;
}


void boundary_condition(float &x, float &y, float &z, float &xmin, float &xmax, float &ymin, float &ymax, float &zmin, float &zmax){
    if (x < xmin)     {x = x + xmax - xmin;}
    else if (x > xmax){x = x - (xmax - xmin);}
    if (y < ymin)     {y = y + ymax - ymin;}
    else if (y > ymax){y = y - (ymax - ymin);}
    if (z < zmin)     {z = z + zmax - zmin;}
    else if (z > zmax){z = z - (zmax - zmin);}
}


void move_cartesian(BootesArray<float> &grain_list, float &min_stoppingtime, float &min_dt, float &max_dt, float &t2nextstop, float &t_frame,
                    CartesianCoord &grid){
    /**
    0  1  2  3   4   5   6  7    8             9
    x, y, z, vx, vy, vz, s, rho, particletype, grainID **/
    #pragma omp parallel for schedule(static)
    for (int grainIND = 0; grainIND < grain_list.shape()[0]; grainIND ++){
        float grain_time_now = 0;           // defined to be time since the start of the loop
        float t_frame_here = t_frame;       // defined to be time since last hydro frame (for interpolation)
        while (grain_time_now < t2nextstop){
            float w[3][3][3];
            int i, j, k;
            TSC(grain_list(grainIND, 0), grain_list(grainIND, 1), grain_list(grainIND, 2),
                grid.xlist, grid.dxlist, grid.ylist, grid.dylist, grid.zlist, grid.dzlist,
                w, i, j, k);
            float gas_vx  = apply_TSC_slope_periodicBC(w, i, j, k, grid.vx,    grid.nx, grid.ny, grid.nz, t_frame_here);
            float gas_vy  = apply_TSC_slope_periodicBC(w, i, j, k, grid.vy,    grid.nx, grid.ny, grid.nz, t_frame_here);
            float gas_vz  = apply_TSC_slope_periodicBC(w, i, j, k, grid.vz,    grid.nx, grid.ny, grid.nz, t_frame_here);
            float grav_x  = apply_TSC_slope_periodicBC(w, i, j, k, grid.gravx, grid.nx, grid.ny, grid.nz, t_frame_here);
            float grav_y  = apply_TSC_slope_periodicBC(w, i, j, k, grid.gravy, grid.nx, grid.ny, grid.nz, t_frame_here);
            float grav_z  = apply_TSC_slope_periodicBC(w, i, j, k, grid.gravz, grid.nx, grid.ny, grid.nz, t_frame_here);
            float gas_rho = apply_TSC_slope_periodicBC(w, i, j, k, grid.rho,   grid.nx, grid.ny, grid.nz, t_frame_here);
            float gas_T   = apply_TSC_slope_periodicBC(w, i, j, k, grid.temp,  grid.nx, grid.ny, grid.nz, t_frame_here);

            float dvx = gas_vx - grain_list(grainIND, 3);
            float dvy = gas_vy - grain_list(grainIND, 4);
            float dvz = gas_vz - grain_list(grainIND, 5);

            float ts = max(min_stoppingtime, grain_list(grainIND, 6) * grain_list(grainIND, 7) / (gas_rho * vth(gas_T)));

            float ax_drag = 1.0 / ts * dvx + grav_x;
            float ay_drag = 1.0 / ts * dvy + grav_y;
            float az_drag = 1.0 / ts * dvz + grav_z;

            float dt_adj_with_ts = max(min_dt, min(max_dt, ts));
            float dt_adj = min(dt_adj_with_ts, t2nextstop - grain_time_now);

            grain_list(grainIND, 0) += grain_list(grainIND, 3) * dt_adj + 0.5 * ax_drag * pow(dt_adj, 2.0);
            grain_list(grainIND, 1) += grain_list(grainIND, 4) * dt_adj + 0.5 * ay_drag * pow(dt_adj, 2.0);
            grain_list(grainIND, 2) += grain_list(grainIND, 5) * dt_adj + 0.5 * az_drag * pow(dt_adj, 2.0);
            grain_list(grainIND, 3) += ax_drag * dt_adj;
            grain_list(grainIND, 4) += ay_drag * dt_adj;
            grain_list(grainIND, 5) += az_drag * dt_adj;

            grain_time_now += dt_adj;
            t_frame_here += dt_adj;
            boundary_condition(grain_list(grainIND, 0), grain_list(grainIND, 1), grain_list(grainIND, 2),
                               grid.xmin, grid.xmax, grid.ymin, grid.ymax, grid.zmin, grid.zmax);
        }
    }
    #pragma omp barrier
}


void move_cartesian_unitless(BootesArray<float> &grain_list, float &min_stoppingtime, float &min_dt, float &max_dt, float &t2nextstop, float &t_frame,
                    CartesianCoord &grid){
    /**
    0  1  2  3   4   5   6  7    8             9
    x, y, z, vx, vy, vz, s, rho, particletype, grainID

    Note: gas_T is gas thermal speed vth = sqrt(8kT/(pi*mu*mH))
    **/
    #pragma omp parallel for schedule(dynamic, 100)
    for (int grainIND = 0; grainIND < grain_list.shape()[0]; grainIND ++){
        //cout << grainIND << endl << flush;
        float grain_time_now = 0;           // defined to be time since the start of the loop
        float t_frame_here = t_frame;       // defined to be time since last hydro frame (for interpolation)
        while (grain_time_now < t2nextstop){
            float w[3][3][3];
            int i, j, k;
            TSC(grain_list(grainIND, 0), grain_list(grainIND, 1), grain_list(grainIND, 2),
                grid.xlist, grid.dxlist, grid.ylist, grid.dylist, grid.zlist, grid.dzlist,
                w, i, j, k);
            float gas_vx  = apply_TSC_slope_periodicBC(w, i, j, k, grid.vx,    grid.nx, grid.ny, grid.nz, t_frame_here);
            float gas_vy  = apply_TSC_slope_periodicBC(w, i, j, k, grid.vy,    grid.nx, grid.ny, grid.nz, t_frame_here);
            float gas_vz  = apply_TSC_slope_periodicBC(w, i, j, k, grid.vz,    grid.nx, grid.ny, grid.nz, t_frame_here);
            float grav_x  = apply_TSC_slope_periodicBC(w, i, j, k, grid.gravx, grid.nx, grid.ny, grid.nz, t_frame_here);
            float grav_y  = apply_TSC_slope_periodicBC(w, i, j, k, grid.gravy, grid.nx, grid.ny, grid.nz, t_frame_here);
            float grav_z  = apply_TSC_slope_periodicBC(w, i, j, k, grid.gravz, grid.nx, grid.ny, grid.nz, t_frame_here);
            float gas_rho = apply_TSC_slope_periodicBC(w, i, j, k, grid.rho,   grid.nx, grid.ny, grid.nz, t_frame_here);
            float gas_T   = apply_TSC_slope_periodicBC(w, i, j, k, grid.temp,  grid.nx, grid.ny, grid.nz, t_frame_here);     // vth thermal speed

            float dvx = gas_vx - grain_list(grainIND, 3);
            float dvy = gas_vy - grain_list(grainIND, 4);
            float dvz = gas_vz - grain_list(grainIND, 5);

            float ts = max(min_stoppingtime, grain_list(grainIND, 6) * grain_list(grainIND, 7) / (gas_rho * gas_T));

            float ax_drag = 1.0 / ts * dvx + grav_x;
            float ay_drag = 1.0 / ts * dvy + grav_y;
            float az_drag = 1.0 / ts * dvz + grav_z;

            float dt_adj_with_ts = max(min_dt, min(max_dt, ts));
            float dt_adj = min(dt_adj_with_ts, t2nextstop - grain_time_now);

            grain_list(grainIND, 0) += grain_list(grainIND, 3) * dt_adj + 0.5 * ax_drag * pow(dt_adj, 2.0);
            grain_list(grainIND, 1) += grain_list(grainIND, 4) * dt_adj + 0.5 * ay_drag * pow(dt_adj, 2.0);
            grain_list(grainIND, 2) += grain_list(grainIND, 5) * dt_adj + 0.5 * az_drag * pow(dt_adj, 2.0);
            grain_list(grainIND, 3) += ax_drag * dt_adj;
            grain_list(grainIND, 4) += ay_drag * dt_adj;
            grain_list(grainIND, 5) += az_drag * dt_adj;

            grain_time_now += dt_adj;
            t_frame_here += dt_adj;
            boundary_condition(grain_list(grainIND, 0), grain_list(grainIND, 1), grain_list(grainIND, 2),
                               grid.xmin, grid.xmax, grid.ymin, grid.ymax, grid.zmin, grid.zmax);
        }
    }
    //#pragma omp barrier
}



#endif

// int main(){ return 0; }

