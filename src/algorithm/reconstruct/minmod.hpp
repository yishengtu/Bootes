#ifndef MINMOD_RECONSTRUCT_HPP_
#define MINMOD_RECONSTRUCT_HPP_

#include "../BootesArray.hpp"
#include "../eos/momentum.hpp"
#include "../../defs.hpp"


class mesh;


void minmod(double &quanp1, double &quan, double &quanm1, double &dx_axis, double &dt, double &Vui, double &acs, double &BquanL, double &BquanR);


void reconstruct_minmod(mesh &m,
                   BootesArray<double> &valsL,
                   BootesArray<double> &valsR,
                   int &x1excess, int &x2excess, int &x3excess,
                   int &axis,
                   int &IMP,
                   double &dt
                   );


#endif
