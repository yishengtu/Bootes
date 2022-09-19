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
    BootesArray<double> valsL;      // boundary left value
    BootesArray<double> valsR;      // boundary right value
    BootesArray<double> fcons;      // flux of conservative variables
    valsL.NewBootesArray(3, NUMCONS, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    valsR.NewBootesArray(3, NUMCONS, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    fcons.NewBootesArray(NUMCONS, 3, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    # pragma acc enter data copyin(valsL, valsR, fcons)
    calc_flux(m, dt, fcons, valsL, valsR);
    #ifdef ENABLE_VISCOSITY
        apply_viscous_flux(m, dt, fcons, m.nu_vis);
    #endif // ENABLE_VISCOSITY
    #if defined(ENABLE_DUSTFLUID)
    // TODO: the nan values probably comes from the fact that v_dust >> v_gas,
    // so the CFL is not satisfied for dust. Periahps the way to get around this is to invoke
    // adaptive time step, for grains which needs to evolve with more time steps
    BootesArray<double> dvalsL;      // boundary left value
    BootesArray<double> dvalsR;      // boundary right value
    BootesArray<double> fdcons;      // flux of dconservative variables
    dvalsL.NewBootesArray(m.NUMSPECIES, 3, NUMCONS - 1, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    dvalsR.NewBootesArray(m.NUMSPECIES, 3, NUMCONS - 1, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    fdcons.NewBootesArray(m.NUMSPECIES, NUMCONS - 1, 3, m.nx3 + 1, m.nx2 + 1, m.nx1 + 1);
    // TODO
    # pragma acc enter data copyin(dvalsL, dvalsR, fdcons)
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
        # pragma acc enter data copyin(stoppingtime_mesh)
        calc_stoppingtimemesh(m, stoppingtime_mesh);

        advect_cons_dust(m, dt, m.NUMSPECIES, fdcons, dvalsL, dvalsR, stoppingtime_mesh);
        #ifdef ENABLE_DUST_GRAINGROWTH
            grain_growth(m, stoppingtime_mesh, dt);
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


