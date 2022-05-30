#ifndef GAS_DRAG_ON_DUST_HPP_
#define GAS_DRAG_ON_DUST_HPP_
#include "../BootesArray.hpp"

class mesh;


double stoppingtime(double &rhodmsize, double &rho, double &pres, double &vth_coeff);


double calc_stoppingtimemesh(mesh &m, BootesArray<double> &stoppingtimemesh);


void apply_gas_drag_on_dust(mesh &m, double &dt, BootesArray<double> &stoppingtimemesh);


#endif // GAS_DRAG_ON_DUST_HPP_

