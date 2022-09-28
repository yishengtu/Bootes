#ifndef COAGULATION_HPP_
#define COAGULATION_HPP_
#include "../../BootesArray.hpp"
#include "../../physical_constants.hpp"


class mesh;


void grain_growth(mesh &m, double dt);


#endif // COAGULATION_HPP_
