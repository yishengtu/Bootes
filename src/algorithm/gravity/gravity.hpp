#ifndef GRAVITY_HPP_
#define GRAVITY_HPP_

#include "../mesh.hpp"
#include "../BootesArray.hpp"
#include "../physical_constants.hpp"
#include "../../defs.hpp"

class gravity{
    public:
        gravity(mesh &m);
        mesh* mpt;
        void pointsource_grav(double m_source, double x1_s, double x2_s, double x3_s);
};


gravity::gravity(mesh &m){
    mpt = &m;
    mpt->Phi_grav.NewBootesArray(m.x3v.shape()[0], m.x2v.shape()[0], m.x1v.shape()[0]);
}


void gravity::pointsource_grav(double m_source, double x1_s, double x2_s, double x3_s){
    for (int kk = mpt->x3s; kk < mpt->x3l; kk ++){
        for (int jj = mpt->x2s; jj < mpt->x2l; jj ++){
            for (int ii = mpt->x1s; ii < mpt->x1l; ii ++){
                double x1_h = mpt->x1v(ii);
                double x2_h = mpt->x2v(jj);
                double x3_h = mpt->x3v(kk);
                double r;
                #if defined(CARTESIAN_COORD)
                    r = sqrt(pow(x1h-x1s, 2) + pow(x2h-x2s, 2) + pow(x3h-x3s, 2));
                #elif defined(SPHERICAL_POLAR_COORD)
                    r = sqrt(pow(x1_h, 2) + pow(x1_s, 2) - 2 * x1_h * x1_s * (sin(x2_h) * sin(x2_s) * cos(x3_h - x3_s) + cos(x2_h) * cos(x2_s)));
                #endif // defined
                mpt->Phi_grav(kk, jj, ii) = G * m_source / r;
            }
        }
    }
}


#endif // GRAVITY_HPP_

