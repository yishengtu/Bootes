#ifndef HYDROGRAV_HPP_
#define HYDROGRAV_HPP_

#include "../../BootesArray.hpp"
#include "../../../defs.hpp"
#include "../../index_def.hpp"


class mesh;


#if defined (ENABLE_GRAVITY)
void apply_grav_source_terms(mesh &m, double &dt);
#endif // defined (enable gravity)


#endif // HYDROGRAV_HPP_
