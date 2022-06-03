#include "timeintegration.hpp"
#include "adv_hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/minmod.hpp"
//#include "reconstruct/const_recon.hpp"
#include "../time_step/time_step.hpp"
#include "../BootesArray.hpp"
#include "../util/util.hpp"
#include "../boundary_condition/apply_bc.hpp"
#include "../index_def.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterm/hydrograv.hpp"

#ifdef ENABLE_DUSTFLUID
    #include "adv_dust.hpp"
    #include "../dust/gas_drag_on_dust.hpp"
    #include "../eos/eos_dust.hpp"
    // #include "../dust/srcterm/dustsrc_term.hpp"
    #include "../boundary_condition/dust/apply_bc_dust.hpp"
#endif // ENABLE_DUSTFLUID


void first_order(mesh &m, double &dt){
    // First order integration
    // (axis, z, y, x)
    /** Step 1: calculate flux **/
    BootesArray<double> valsL;      // boundary left value
    BootesArray<double> valsR;      // boundary right value
    BootesArray<double> fcons;      // flux of conservative variables
    valsL.NewBootesArray(3, NUMCONS, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    valsR.NewBootesArray(3, NUMCONS, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    fcons.NewBootesArray(NUMCONS, 3, m.cons.shape()[1] + 1, m.cons.shape()[2] + 1, m.cons.shape()[3] + 1);
    calc_flux(m, dt, fcons, valsL, valsR);
    #if defined(ENABLE_DUSTFLUID)
    // TODO: the nan values probably comes from the fact that v_dust >> v_gas,
    // so the CFL is not satisfied for dust. Periahps the way to get around this is to invoke
    // adaptive time step, for grains which needs to evolve with more time steps
    BootesArray<double> dvalsL;      // boundary left value
    BootesArray<double> dvalsR;      // boundary right value
    BootesArray<double> fdcons;      // flux of dconservative variables
    dvalsL.NewBootesArray(m.NUMSPECIES, 3, NUMCONS - 1, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    dvalsR.NewBootesArray(m.NUMSPECIES, 3, NUMCONS - 1, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    fdcons.NewBootesArray(m.NUMSPECIES, NUMCONS - 1, 3, m.dcons.shape()[2] + 1, m.dcons.shape()[3] + 1, m.dcons.shape()[4] + 1);
    calc_flux_dust(m, dt, m.NUMSPECIES, fdcons, dvalsL, dvalsR);
    #endif

    /** step 2: hydro: time integrate to update CONSERVATIVE variables, solve Riemann Problem **/
    /** step 2.1: hydro **/
    advect_cons(m, dt, fcons, valsL, valsR);
    /** step 2.2: hydro source **/
    #if defined (ENABLE_GRAVITY)
        apply_grav_source_terms(m, dt);
    #endif
    /** step 3: dust: time integrate to update CONSERVATIVE variables, solve Riemann Problem **/
    /** step 3.1: dust **/
    #ifdef ENABLE_DUSTFLUID
        BootesArray<double> stoppingtime_mesh;
        stoppingtime_mesh.NewBootesArray(m.NUMSPECIES, m.x3v.shape()[0], m.x2v.shape()[0], m.x1v.shape()[0]);
        double min_stoptime = calc_stoppingtimemesh(m, stoppingtime_mesh);

        advect_cons_dust(m, dt, m.NUMSPECIES, fdcons, dvalsL, dvalsR, stoppingtime_mesh);
        /** step 3.2: dust source **/
        /* all source terms are taken care in one function,
           because if ts too small the values are directly replaced with terminal velocity approximation.
        */
        // apply_source_terms_dust(m, dt, stoppingtime_mesh);
    #endif // ENABLE_DUSTFLUID

    /** step 4: protections **/
    #if defined (PROTECTION_PROTECTION)
        protection(m);
    #endif // defined(PROTECTION_PROTECTION)
    #ifdef DUST_PROTECTION
        protection_dust(m);
    #endif // DUST_PROTECTION

    /** step 5: use E.O.S. and relations to get primitive variables. **/
    cons_to_prim(m);
    #ifdef ENABLE_DUSTFLUID
    cons_to_prim_dust(m);
    #endif // ENABLE_DUSTFLUID

    /** step 6: apply boundary conditions **/
    apply_boundary_condition(m);
    #ifdef ENABLE_DUSTFLUID
    apply_boundary_condition_dust(m);
    #endif
}


