#ifndef ADV_HYDRO_HPP_
#define ADV_HYDRO_HPP_

#include "../BootesArray.hpp"
#include "../../defs.hpp"
#include "../index_def.hpp"


class mesh;


void calc_flux(mesh &m, double dt);


void advect_cons(mesh &m, double dt);


#ifdef DENSITY_PROTECTION
void protection(mesh &m);
#endif // DENSITY_PROTECTION

#ifdef ENABLE_TEMPERATURE_PROTECTION
void temperature_protection(mesh &m);
#endif // ENABLE_TEMPERATURE_PROTECTION

#endif // ADV_HYDRO_HPP_
