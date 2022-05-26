#include "standard_bc.hpp"
#include "periodic_bc.hpp"
#include "reflective_bc.hpp"
#include "spherical_polar_pole.hpp"
#include "../mesh/mesh.hpp"


void apply_boundary_condition(mesh &m){
    standard_boundary_condition_x1i(m.cons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    standard_boundary_condition_x1i(m.prim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    standard_boundary_condition_x1o(m.cons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    standard_boundary_condition_x1o(m.prim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    standard_boundary_condition_x2i(m.cons, m.x1s, m.x1l, m.ng1,
                                        m.x2s, m.x2l, m.ng2,
                                        m.x3s, m.x3l, m.ng3);
    standard_boundary_condition_x2i(m.prim, m.x1s, m.x1l, m.ng1,
                                        m.x2s, m.x2l, m.ng2,
                                        m.x3s, m.x3l, m.ng3);

    standard_boundary_condition_x2o(m.cons, m.x1s, m.x1l, m.ng1,
                                        m.x2s, m.x2l, m.ng2,
                                        m.x3s, m.x3l, m.ng3);
    standard_boundary_condition_x2o(m.prim, m.x1s, m.x1l, m.ng1,
                                        m.x2s, m.x2l, m.ng2,
                                        m.x3s, m.x3l, m.ng3);

    standard_boundary_condition_x3i(m.cons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    standard_boundary_condition_x3i(m.prim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    standard_boundary_condition_x3o(m.prim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    standard_boundary_condition_x3o(m.cons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
}

