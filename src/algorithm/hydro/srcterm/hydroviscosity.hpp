#ifndef HYDROVISCOSITY_HPP_
#define HYDROVISCOSITY_HPP_

#include "../../BootesArray.hpp"
#include "../../../defs.hpp"
#include "../../index_def.hpp"


class mesh;


#if defined (ENABLE_VISCOSITY)
void apply_viscous_flux(mesh &m, double &dt, BootesArray<double> &fcons, BootesArray<double> &nu_vis);
#endif // defined (ENABLE_VISCOSITY)


#endif // HYDROVISCOSITY_HPP_
