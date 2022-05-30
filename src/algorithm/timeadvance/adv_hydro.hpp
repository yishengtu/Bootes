#ifndef ADV_HYDRO_HPP_
#define ADV_HYDRO_HPP_

#include "../BootesArray.hpp"
#include "../../defs.hpp"
#include "../index_def.hpp"


class mesh;


void calc_flux(mesh &m, double &dt, BootesArray<double> &fcons, BootesArray<double> &valsL, BootesArray<double> &valsR);


void advect_cons(mesh &m, double &dt, BootesArray<double> &fcons, BootesArray<double> &valsL, BootesArray<double> &valsR);


#ifdef PROTECTION_PROTECTION
void protection(mesh &m);
#endif // PROTECTION_PROTECTION

#if defined (ENABLE_GRAVITY)
void apply_grav_source_terms(mesh &m, double &dt);
#endif // defined (enable gravity)

void first_order(mesh &m, double &dt);


#endif // ADV_HYDRO_HPP_
