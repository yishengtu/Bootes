#ifndef ADV_DUST_HPP_
#define ADV_DUST_HPP_

#include "../BootesArray.hpp"
#include "../index_def.hpp"
#include "../../defs.hpp"

class mesh;


void calc_flux_dust(mesh &m, double &dt, int &NUMSPECIES, BootesArray<double> &fdcons, BootesArray<double> &valsL, BootesArray<double> &valsR);


void advect_cons_dust(mesh &m, double &dt, int &NUMSPECIES, BootesArray<double> &fdcons, BootesArray<double> &valsL, BootesArray<double> &valsR, BootesArray<double> &stoppingtimemesh);


#ifdef ENABLE_GRAVITY
void apply_grav_source_terms_dust(mesh &m, double &dt);
#endif // ENABLE_GRAVITY

void first_order_dust(mesh &m, double &dt);

#endif // ADV_DUST_HPP_
