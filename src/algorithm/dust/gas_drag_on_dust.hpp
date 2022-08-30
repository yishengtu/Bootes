#ifndef GAS_DRAG_ON_DUST_HPP_
#define GAS_DRAG_ON_DUST_HPP_
#include "../BootesArray.hpp"

class mesh;


double stoppingtime(double &rhodmsize, double &rho, double &pres, double &vth_coeff);


void calc_stoppingtimemesh(mesh &m, BootesArray<double> &stoppingtimemesh);


double find_smallest_stoppingtime(mesh &m, BootesArray<double> &stoppingtimemesh);


#endif // GAS_DRAG_ON_DUST_HPP_

