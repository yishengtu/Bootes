#ifndef STANDARD_BC_HPP_
#define STANDARD_BC_HPP_
#include "../BootesArray.hpp"
#include "../index_def.hpp"


void standard_boundary_condition_x1i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3);

void standard_boundary_condition_x1o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3);

void standard_boundary_condition_x2i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3);

void standard_boundary_condition_x2o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3);

void standard_boundary_condition_x3i(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3);


void standard_boundary_condition_x3o(BootesArray<double> &quan, int &x1s, int &x1l, int &ng1,
                                                            int &x2s, int &x2l, int &ng2,
                                                            int &x3s, int &x3l, int &ng3);

#endif // STANDARD_BC_HPP_

