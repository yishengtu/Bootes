#ifndef EOS_DUST_HPP_
#define EOS_DUST_HPP_
#include "../../defs.hpp"

#include <cmath>


class mesh;


#ifdef ENABLE_DUSTFLUID
void cons_to_prim_dust(mesh &m);
#endif // ENABLE_DUSTFLUID


#endif // EOS_DUST_HPP_
