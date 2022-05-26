#ifndef GRAVITY_HPP_
#define GRAVITY_HPP_

#include "../BootesArray.hpp"
#include "../physical_constants.hpp"
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

        void setup_Phimesh(int &tot_nx3, int &tot_nx2, int &tot_nx1);
        void pointsource_grav(mesh &m, double m_source, double x1_s, double x2_s, double x3_s);
        void calc_surface_vals(mesh &m);

};

#endif // GRAVITY_HPP_

