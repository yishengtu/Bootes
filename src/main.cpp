#include "algorithm/mesh/mesh.hpp"
#include "algorithm/util/util.hpp"
#include "algorithm/time_step/time_step.hpp"
#include "algorithm/hydro/hllc.hpp"
#include "algorithm/timeadvance/timeintegration.hpp"
#include "algorithm/hydro/donercell.hpp"
#include "algorithm/inoutput/output.hpp"
#include "algorithm/inoutput/input.hpp"
#include "algorithm/index_def.hpp"
#include "algorithm/boundary_condition/apply_bc.hpp"
#include "algorithm/eos/eos.hpp"
#include "algorithm/physical_constants.hpp"
#include "defs.hpp"
#include <chrono>
#include <omp.h>

#if defined(ENABLE_GRAVITY)
    #include "algorithm/gravity/gravity.hpp"
#endif // defined(ENABLE_GRAVITY)

#ifdef ENABLE_DUSTFLUID
    #include "algorithm/eos/eos_dust.hpp"
    #include "algorithm/boundary_condition/dust/apply_bc_dust.hpp"
#endif // ENABLE_DUSTFLUID

#ifdef DEBUG
    #include "algorithm/util/checkok.hpp"
#endif // DEBUG

#include "setup/test_disksim.cpp"

void doloop(double &ot, double &next_exit_loop_time, mesh &m, double &CFL){
    int loop_cycle = 0;
    while (ot < next_exit_loop_time){
        double dt = timestep(m, CFL);
        dt = min(dt, next_exit_loop_time - ot);
        cout << "\t integrate cycle: " << loop_cycle << "\t time: " << ot << "\t dt: " << dt << endl << flush;
        // step 1: evolve the hydro by dt
        first_order(m, dt);
        // step 2: update other fields
        // step 2.1: calculate source terms
        // step 2.1.1: gravity
        //#if defined (ENABLE_GRAVITY)
        //    m.grav->pointsource_grav(m, 1.e4, 0, 0, 0);
        //    m.grav->calc_surface_vals(m);
        //#endif // defined (ENABLE_GRAVITY)

        // step 3: work after loop
        work_after_loop(m, dt);

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

        #ifdef DEBUG
            /** check the values are fine in program **/
            #ifdef DEBUG
            int stat = check_ok(m);
            #endif // DEBUG
            if (stat == 1){
                throw 1;
            }
        #endif // DEBUG

        // last step: iterate counter
        ot += dt;
        loop_cycle += 1;
    }
}


int main(){
    /** start timer **/
    auto start = std::chrono::steady_clock::now();

    /** >>> these lines exist for output use only **/
    const int outputIDN = static_cast<int>(IDN);
    const int outputIM1 = static_cast<int>(IM1);
    const int outputIM2 = static_cast<int>(IM2);
    const int outputIM3 = static_cast<int>(IM3);
    const int outputIEN = static_cast<int>(IEN);
    const int outputIDP = static_cast<int>(IDP);
    const int outputIV1 = static_cast<int>(IV1);
    const int outputIV2 = static_cast<int>(IV2);
    const int outputIV3 = static_cast<int>(IV3);
    const int outputIPN = static_cast<int>(IPN);
    /** << end for output only << **/

    /** read in necessary information from input file **/
    input_file finput("input.txt");

    double CFL = finput.getDouble("CFL");

    /** setup grid and initial condition **/
    mesh m;
    m.pconst.setup_physical_constants(finput);
    double gamma_hydro = finput.getDouble("gamma_hydro");

    int dim       = finput.getInt("dimension");
    double x1min  = finput.getDouble("x1min");
    double x1max  = finput.getDouble("x1max");
    int nx1       = finput.getInt("nx1");
    double x2min  = finput.getDouble("x2min");
    double x2max  = finput.getDouble("x2max");
    int nx2       = finput.getInt("nx2");
    double x3min  = finput.getDouble("x3min");
    double x3max  = finput.getDouble("x3max");
    int nx3       = finput.getInt("nx3");
    double ratio1 = pow(x1max / x1min, (1./ (double) nx1));
    // determine number of ghost zones
    int ng1, ng2, ng3;
    if      (dim == 1){ ng1 = 2; ng2 = 0; ng3 = 0; }
    else if (dim == 2){ ng1 = 2; ng2 = 2; ng3 = 0; }
    else if (dim == 3){ ng1 = 2; ng2 = 2; ng3 = 2; }
    else { cout << "dimension not recognized! " << endl << flush; throw 1; }
    #if defined(CARTESIAN_COORD)
    m.SetupCartesian(dim,
                      x1min, x1max, nx1, ng1,                       // ax1
                      x2min, x2max, nx2, ng2,                       // ax2
                      x3min, x3max, nx3, ng3                        // ax3
                      );
    #elif defined(SPHERICAL_POLAR_COORD)
    m.SetupSphericalPolar(dim,
                          x1min, x1max, nx1, ratio1, ng1,                       // ax1
                          x2min, x2max, nx2,         ng2,                       // ax2
                          x3min, x3max, nx3,         ng3                        // ax3
                          );
    #endif // defined (COORDINATE)
    m.hydro_gamma = gamma_hydro;
    m.vth_coeff = 8.0 / M_PI * gamma_hydro;         // for calculating gas thermal speed

    #ifdef ENABLE_DUSTFLUID
        setup_dust(m, finput);          // fill in m.GrainEdgeList, m.GrainSizeList and m.NUMSPECIES
        m.setupDustFluidMesh(m.NUMSPECIES);
        for (int ii = 0; ii < m.NUMSPECIES; ii ++){
            cout << m.GrainSizeList(ii) << endl << flush;
        }
        m.GrainSizeTimesGrainDensity.NewBootesArray(m.NUMSPECIES);
        m.GrainMassList.NewBootesArray(m.NUMSPECIES);
        for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
            m.GrainSizeTimesGrainDensity(specIND) = m.GrainSizeList(specIND) * m.rhodm;
            m.GrainMassList(specIND) = 4. / 3. * M_PI * pow(m.GrainSizeList(specIND), 3) * m.rhodm;
        }
    #endif // ENABLE_DUSTFLUID
    /** setup initial condition **/
    setup(m, finput);   // setup according to the input file

    cons_to_prim(m);
    apply_boundary_condition(m);
    #ifdef ENABLE_DUSTFLUID
    cons_to_prim_dust(m);
    apply_boundary_condition_dust(m);
    #endif
    /** initialize simulation parameters **/
    double t_tot = finput.getDouble("t_tot");
    double output_dt = finput.getDouble("output_dt");
    string foutput_root = finput.getString("foutput_root");
    string foutput_pre  = finput.getString("foutput_pre");
    string foutput_aft  = finput.getString("foutput_aft");

    /** initialize time and cycle trackings **/
    int frame = 0;
    double ot = 0;
    int cycle = 0;

    /** initialize decisions in loop **/
    double next_output_time = ot;
    double next_exit_loop_time = next_output_time;
    bool det_output = false;
    bool det_doloop = false;

    std::cout << "setup complete" << std::endl << flush;
    /** main loop **/
    while (ot < t_tot){
        // step 1: determine when to exit the time integration loop
        next_exit_loop_time = min(next_output_time, t_tot);
        // step 2: determine what needs to be done
        if (ot >= next_output_time){
            det_output = true;
            det_doloop = false;
        }
        else{
            det_output = false;
            det_doloop = true;
        }
        // step 3: do what needs to be done
        if (det_doloop){
            doloop(ot, next_exit_loop_time, m, CFL);
        }
        if (det_output){
            Output output(foutput_root + foutput_pre + "." + choosenumber(frame) + "." + foutput_aft, 'w');
            output.writeattribute<double>(&ot, "time", H5::PredType::NATIVE_DOUBLE, 1);
            output.writeattribute<int>(&frame, "frame", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<int>(&cycle, "main_cycle", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<int>(&m.x1s, "x1s", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<int>(&m.x1l, "x1l", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<int>(&m.x2s, "x2s", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<int>(&m.x2l, "x2l", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<int>(&m.x3s, "x3s", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<int>(&m.x3l, "x3l", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIDN, "rhoIND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIM1, "mo1IND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIM2, "mo2IND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIM3, "mo3IND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIEN, "eneIND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIDP, "denIND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIV1, "ve1IND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIV2, "ve2IND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIV3, "ve3IND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<const int>(&outputIPN, "prsIND", H5::PredType::NATIVE_INT32, 1);
            output.writeattribute<double>(&m.hydro_gamma, "hydro_gamma", H5::PredType::NATIVE_DOUBLE, 1);
            output.write1Ddataset(m.x1v, "x1v", H5::PredType::NATIVE_DOUBLE);
            output.write1Ddataset(m.x2v, "x2v", H5::PredType::NATIVE_DOUBLE);
            output.write1Ddataset(m.x3v, "x3v", H5::PredType::NATIVE_DOUBLE);
            output.write1Ddataset(m.x1f, "x1f", H5::PredType::NATIVE_DOUBLE);
            output.write1Ddataset(m.x2f, "x2f", H5::PredType::NATIVE_DOUBLE);
            output.write1Ddataset(m.x3f, "x3f", H5::PredType::NATIVE_DOUBLE);
            output.write3Ddataset(m.dx1p, "dx1p", H5::PredType::NATIVE_DOUBLE);
            output.write3Ddataset(m.dx2p, "dx2p", H5::PredType::NATIVE_DOUBLE);
            output.write3Ddataset(m.dx3p, "dx3p", H5::PredType::NATIVE_DOUBLE);
            output.write4Ddataset(m.prim, "prim", H5::PredType::NATIVE_DOUBLE);
            output.write4Ddataset(m.cons, "cons", H5::PredType::NATIVE_DOUBLE);
            #if defined(ENABLE_GRAVITY)
            output.write3Ddataset(m.grav->Phi_grav, "Phi", H5::PredType::NATIVE_DOUBLE);
            output.write3Ddataset(m.grav->Phi_grav_x1surface, "Phi_x1s", H5::PredType::NATIVE_DOUBLE);
            output.write3Ddataset(m.grav->Phi_grav_x2surface, "Phi_x2s", H5::PredType::NATIVE_DOUBLE);
            output.write3Ddataset(m.grav->Phi_grav_x3surface, "Phi_x3s", H5::PredType::NATIVE_DOUBLE);
            #endif
            #if defined(ENABLE_DUSTFLUID)
            output.write5Ddataset(m.dcons, "dcons", H5::PredType::NATIVE_DOUBLE);
            output.write5Ddataset(m.dprim, "dprim", H5::PredType::NATIVE_DOUBLE);
            #endif
            if (m.UserScalers.checkallocated()){
            output.write1Ddataset(m.UserScalers, "UserScalers", H5::PredType::NATIVE_DOUBLE);
            }
            output.close();
            double elasped = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() / 1000.;
            std::cout << "Output frame " << frame << '\t' << "Elapsed real time =" << elasped << " seconds" << std::endl;
            frame += 1;
            next_output_time += output_dt;
        }
        cycle += 1;
        cout << "main cycle: " << cycle << "    time: " << ot << endl << flush;
    }
    return 0;
}
