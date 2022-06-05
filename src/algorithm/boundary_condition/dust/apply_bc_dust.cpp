#include "apply_bc_dust.hpp"
#include "standard_bc_dust.hpp"
//#include "periodic_bc_dust.hpp"
#include "reflective_bc_dust.hpp"
//#include "spherical_polar_pole_dust.hpp"
#include "../../mesh/mesh.hpp"


void apply_boundary_condition_dust(mesh &m){
    dust_standard_boundary_condition_x1i(m.dcons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    dust_standard_boundary_condition_x1i(m.dprim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    dust_standard_boundary_condition_x1o(m.dcons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    dust_standard_boundary_condition_x1o(m.dprim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    dust_reflective_boundary_condition_x2i(m.dcons, m.x1s, m.x1l, m.ng1,
                                        m.x2s, m.x2l, m.ng2,
                                        m.x3s, m.x3l, m.ng3);
    dust_reflective_boundary_condition_x2i(m.dprim, m.x1s, m.x1l, m.ng1,
                                        m.x2s, m.x2l, m.ng2,
                                        m.x3s, m.x3l, m.ng3);

    dust_reflective_boundary_condition_x2o(m.dcons, m.x1s, m.x1l, m.ng1,
                                        m.x2s, m.x2l, m.ng2,
                                        m.x3s, m.x3l, m.ng3);
    dust_reflective_boundary_condition_x2o(m.dprim, m.x1s, m.x1l, m.ng1,
                                        m.x2s, m.x2l, m.ng2,
                                        m.x3s, m.x3l, m.ng3);

    dust_standard_boundary_condition_x3i(m.dcons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    dust_standard_boundary_condition_x3i(m.dprim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);

    dust_standard_boundary_condition_x3o(m.dcons, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
    dust_standard_boundary_condition_x3o(m.dprim, m.x1s, m.x1l, m.ng1,
                                           m.x2s, m.x2l, m.ng2,
                                           m.x3s, m.x3l, m.ng3);
}

