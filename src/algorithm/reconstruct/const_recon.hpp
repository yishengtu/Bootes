#ifndef CONST_RECONSTRUCT_HPP_
#define CONST_RECONSTRUCT_HPP_

#include "../BootesArray.hpp"
#include "../eos/momentum.hpp"
#include "../../defs.hpp"


class mesh;


void reconstruct_const(mesh &m,
                       BootesArray<double> &valsL,
                       BootesArray<double> &valsR,
                       int &x1excess, int &x2excess, int &x3excess,
                       int &axis,
                       int &IMP,
                       double &dt
                       );


#endif
