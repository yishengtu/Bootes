#include "timeintegration.hpp"
#include "adv_hydro.hpp"
#include "../mesh/mesh.hpp"
//#include "../reconstruct/minmod.hpp"
//#include "../reconstruct/MUSCL_Hancock.hpp"
//#include "reconstruct/const_recon.hpp"
#include "../time_step/time_step.hpp"
#include "../BootesArray.hpp"
#include "../util/util.hpp"
#include "../boundary_condition/apply_bc.hpp"
#include "../index_def.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterm/hydrograv.hpp"

#ifdef ENABLE_VISCOSITY
    #include "../hydro/srcterm/hydroviscosity.hpp"
#endif // ENABLE_VISCOSITY

#ifdef ENABLE_DUSTFLUID
    #include "adv_dust.hpp"
    #include "../dust/gas_drag_on_dust.hpp"
    #include "../eos/eos_dust.hpp"
    // #include "../dust/srcterm/dustsrc_term.hpp"
    #include "../boundary_condition/dust/apply_bc_dust.hpp"
    #ifdef ENABLE_DUST_GRAINGROWTH
        #include "../dust/graingrowth/coagulation.hpp"
    #endif // ENABLE_DUST_GRAINGROWTH
#endif // ENABLE_DUSTFLUID


void first_order(mesh &m, double &dt){
    // First order integration
    // (axis, z, y, x)
    /** Step 1: calculate flux **/
    calc_flux(m, dt, m.fcons, m.valsL, m.valsR);
    #ifdef ENABLE_VISCOSITY
        apply_viscous_flux(m, dt, m.fcons, m.nu_vis);
    #endif // ENABLE_VISCOSITY
    #if defined(ENABLE_DUSTFLUID)
    // TODO: the nan values probably comes from the fact that v_dust >> v_gas,
    // so the CFL is not satisfied for dust. Periahps the way to get around this is to invoke
    // adaptive time step, for grains which needs to evolve with more time steps
    calc_flux_dust(m, dt, m.NUMSPECIES, m.fdcons, m.dvalsL, m.dvalsR);
    // Another way to do this is by creating a structured data region by putting {} around a region.
    #endif

    /** step 2: hydro: time integrate to update CONSERVATIVE variables, solve Riemann Problem **/
    /** step 2.1: hydro **/
    advect_cons(m, dt, m.fcons, m.valsL, m.valsR);
    /** step 2.2: hydro source **/
    #if defined (ENABLE_GRAVITY)
        apply_grav_source_terms(m, dt);
    #endif
    /** step 3: dust: time integrate to update CONSERVATIVE variables, solve Riemann Problem **/
    /** step 3.1: dust **/
    #ifdef ENABLE_DUSTFLUID
        calc_stoppingtimemesh(m, m.stoppingtime_mesh);

        advect_cons_dust(m, dt, m.NUMSPECIES, m.fdcons, m.dvalsL, m.dvalsR, m.stoppingtime_mesh);
        #ifdef ENABLE_DUST_GRAINGROWTH
            grain_growth(m, m.stoppingtime_mesh, dt);
        #endif // ENABLE_DUST_GRAINGROWTH
    #endif // ENABLE_DUSTFLUID

    /** step 4: protections **/
    #if defined (DENSITY_PROTECTION)
        protection(m, m.minDensity);
    #endif // defined(DENSITY_PROTECTION)
    #ifdef ENABLE_TEMPERATURE_PROTECTION
        temperature_protection(m, m.minTemp);
    #endif // ENABLE_TEMPERATURE_PROTECTION
}


