#ifndef MUSCL_RECONSTRUCT_HPP_
#define MUSCL_RECONSTRUCT_HPP_

#include "../BootesArray.hpp"


class mesh;


void reconstruct_MHM(mesh &m,
                   BootesArray<double> &valsL,
                   BootesArray<double> &valsR,
                   int &x1excess, int &x2excess, int &x3excess,
                   int &axis,
                   int &IMP,
                   double &dt
                   );


#endif
