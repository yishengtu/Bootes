#ifndef SPHERICAL_POLAR_POLE_BC_HPP_
#define SPHERICAL_POLAR_POLE_BC_HPP_
#include "../BootesArray.hpp"
#include "../index_def.hpp"


void sph_polar_pole_boundary_condition_x2i(BootesArray<double> &quan, int x1s, int x1l, int ng1,
                                                                      int x2s, int x2l, int ng2,
                                                                      int x3s, int x3l, int ng3);

void sph_polar_pole_boundary_condition_x2o(BootesArray<double> &quan, int x1s, int x1l, int ng1,
                                                                      int x2s, int x2l, int ng2,
                                                                      int x3s, int x3l, int ng3);


#endif // STANDARD_BC_HPP_

