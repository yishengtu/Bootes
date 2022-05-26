#ifndef APPLY_BC_HPP_
#define APPLY_BC_HPP_

#include "standard_bc.hpp"
#include "periodic_bc.hpp"
#include "reflective_bc.hpp"
#include "spherical_polar_pole.hpp"


class mesh;


void apply_boundary_condition(mesh &m);


#endif // APPLY_BC_HPP_

