#include "apply_bc.hpp"
#include "../mesh/mesh.hpp"


void apply_boundary_condition(mesh &m){
    outflow_boundary_condition_x1i(m.cons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    outflow_boundary_condition_x1i(m.prim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    outflow_boundary_condition_x1o(m.cons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    outflow_boundary_condition_x1o(m.prim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    sph_polar_pole_boundary_condition_x2i(m.cons, m.x1s, m.x1l, m.ng1,
                                          m.x2s, m.x2l, m.ng2,
                                          m.x3s, m.x3l, m.ng3);
    sph_polar_pole_boundary_condition_x2i(m.prim, m.x1s, m.x1l, m.ng1,
                                          m.x2s, m.x2l, m.ng2,
                                          m.x3s, m.x3l, m.ng3);

    sph_polar_pole_boundary_condition_x2o(m.cons, m.x1s, m.x1l, m.ng1,
                                          m.x2s, m.x2l, m.ng2,
                                          m.x3s, m.x3l, m.ng3);
    sph_polar_pole_boundary_condition_x2o(m.prim, m.x1s, m.x1l, m.ng1,
                                          m.x2s, m.x2l, m.ng2,
                                          m.x3s, m.x3l, m.ng3);

    periodic_boundary_condition_x3i(m.cons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    periodic_boundary_condition_x3i(m.prim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    periodic_boundary_condition_x3o(m.prim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    periodic_boundary_condition_x3o(m.cons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
}

