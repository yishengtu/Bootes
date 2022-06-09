#ifndef GRAVITY_HPP_
#define GRAVITY_HPP_

#include "../BootesArray.hpp"
#include "../../defs.hpp"


class mesh;


class gravity{
    public:
        gravity();

        /** grav **/
        BootesArray<double> Phi_grav;
        BootesArray<double> Phi_grav_x1surface;
        BootesArray<double> Phi_grav_x2surface;
        BootesArray<double> Phi_grav_x3surface;
        BootesArray<double> grav_x1;
        BootesArray<double> grav_x2;
        BootesArray<double> grav_x3;

        void setup_Phimesh(int &tot_nx3, int &tot_nx2, int &tot_nx1);
        void zero_gravity(mesh &m);
        void add_self_grav(mesh &m);
        void add_pointsource_grav(mesh &m, double &m_source, double &x1_s, double &x2_s, double &x3_s);
        void calc_surface_vals(mesh &m);
        void boundary_grav(mesh &m);
        void self_grav(mesh &m);

};

#endif // GRAVITY_HPP_

