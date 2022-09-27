#include "coagulation.hpp"
#include "graingrowthmodel.hpp"
#include "../../util/util.hpp"
#include "../../index_def.hpp"
#include "../../mesh/mesh.hpp"
#include "../../../defs.hpp"
#include "../terminalvel.hpp"
#include <cmath>


#pragma acc routine vector
void grain_growth_one_cell_stick(double num[],
                           double vr[], double vtheta[], double vphi[], double num_here[], double Mmat[],
                           BootesArray<double> &grain_size_list, BootesArray<double> &grain_mass_list, double dt, int NUM_SPECIES){
    double dt_here = dt;
#pragma acc loop vector
    for (int i = 0; i < NUM_SPECIES; i++){num_here[i] = num[i];}
    bool redo = true;
    double dt_tot = 0;
    redo = false;   // set it to false first
#pragma acc loop vector //collapse (2)
    for (int j = 0; j < NUM_SPECIES; ++ j){
        for (int k = j; k < NUM_SPECIES; ++ k){
            double dv_bulk = sqrt(pow((vr[j] - vr[k]), 2) + pow((vtheta[j] - vtheta[k]), 2) + pow((vphi[j] - vphi[k]), 2));
            // double dv_vortex = dv_ormel(grain_size_list(j), grain_size_list(k), rhogas, tempgas);
            double dv = dv_bulk;

            double K = dv * M_PI * pow(grain_size_list(j) + grain_size_list(k), 2.0);

            // replace the for statement here with new method, use V1.8's for if needed
            // gain via coagulation
            double numjtimesnumk = num_here[j] * num_here[k];
            int cog_res = 0;        // the following loop is equivalent to searchsorted in numpy
            while (cog_res < NUM_SPECIES){
                if (grain_mass_list(cog_res) >= (grain_mass_list(j) + grain_mass_list(k))){
                    break;
                }
                cog_res += 1;
            }
            // int cog_res = searchsorted(grain_mass_list(j) + grain_mass_list(k), grain_mass_list) - 1;
            if (j == k){
                if (cog_res == NUM_SPECIES){
#pragma acc atomic update
                    Mmat[cog_res - 1] += 0.5 * (K * grain_mass_list(k) + K * grain_mass_list(j)) / grain_mass_list(cog_res - 1) * numjtimesnumk;
                }
                else{
                    double eps = (grain_mass_list(j) + grain_mass_list(k) - grain_mass_list(cog_res - 1)) / (grain_mass_list(cog_res) - grain_mass_list(cog_res - 1));
#pragma acc atomic update
                    Mmat[cog_res]     += 0.5 * (K * grain_mass_list(k) + K * grain_mass_list(j)) / grain_mass_list(cog_res) * numjtimesnumk * eps;
#pragma acc atomic update
                    Mmat[cog_res - 1] += 0.5 * (K * grain_mass_list(k) + K * grain_mass_list(j)) / grain_mass_list(cog_res - 1) * numjtimesnumk * (1.0 - eps);
                }
            }
            else{
                if (cog_res == NUM_SPECIES){
#pragma acc atomic update
                    Mmat[cog_res - 1] += (K * grain_mass_list(k) + K * grain_mass_list(j)) / grain_mass_list(cog_res - 1)* numjtimesnumk;
                }
                else{
                    double eps =  (grain_mass_list(j) + grain_mass_list(k) - grain_mass_list(cog_res - 1)) / (grain_mass_list(cog_res) - grain_mass_list(cog_res - 1));
#pragma acc atomic update
                    Mmat[cog_res]     += (K * grain_mass_list(k) + K * grain_mass_list(j)) / grain_mass_list(cog_res)* numjtimesnumk * eps;
#pragma acc atomic update
                    Mmat[cog_res - 1] += (K * grain_mass_list(k) + K * grain_mass_list(j)) / grain_mass_list(cog_res - 1)* numjtimesnumk * (1.0 - eps);
                }
            }
            // lost via coagulation
#pragma acc atomic update
            Mmat[k] -= K * numjtimesnumk;
            if (k != j){
#pragma acc atomic update
                Mmat[j] -= K * numjtimesnumk;
            }
            // break if there is a problem
            /*
            if (isnan(Mmat(0))){
            }
            if (isnan(Mmat(j))){
            }
            if (isnan(Mmat(k))){
            }
            */
        }
    }
    dt_tot += dt_here;
#pragma acc loop vector
    for (int i = 0; i < NUM_SPECIES; ++ i){
        num_here[i] += Mmat[i] * dt_here;
        if (num_here[i] < 0)
        {
            num_here[i] = 0;
        }
    }
    //cout << dt_tot << '\t' << min_dt << '\t' << dt << '\n' << flush;
#pragma acc loop vector
    for (int i = 0; i < NUM_SPECIES; i++){num[i] = num_here[i];}
}


#pragma acc routine vector
void grain_growth_one_cell(double num[],
                           double vr[], double vtheta[], double vphi[], double num_here[], double Mmat[],
                           BootesArray<double> &grain_size_list, BootesArray<double> &grain_mass_list, double dt, int NUM_SPECIES){
    double dt_here = dt;
    for (int i = 0; i < NUM_SPECIES; i++){num_here[i] = num[i];}
    bool redo = true;
    double dt_tot = 0;
    while (redo || dt_tot < dt)
    {
        redo = false;   // set it to false first
        // BootesArray<double> Mmat; Mmat.NewBootesArray(NUM_SPECIES);         // initialize with 0s
#pragma acc loop vector
        for (int k = 0; k < NUM_SPECIES; k++){
            Mmat[k] = 0;
        }
#pragma acc loop vector collapse (2)
        for (int j = 0; j < NUM_SPECIES; ++ j){
            for (int k = j; k < NUM_SPECIES; ++ k){
                double dv_bulk = sqrt(pow((vr[j] - vr[k]), 2) + pow((vtheta[j] - vtheta[k]), 2) + pow((vphi[j] - vphi[k]), 2));
                // double dv_vortex = dv_ormel(grain_size_list(j), grain_size_list(k), rhogas, tempgas);
                double dv = dv_bulk;

                double KL1[2];                  //K1 = KL1[0], L1 = KL1[1]
                double KL2[2];
                grain_growth_model_stick(grain_size_list(j), grain_size_list(k), dv, KL1);
                grain_growth_model_stick(grain_size_list(k), grain_size_list(j), dv, KL2);

                // replace the for statement here with new method, use V1.8's for if needed
                // gain via coagulation
                double numjtimesnumk = num_here[j] * num_here[k];
                int cog_res = 0;        // the following loop is equivalent to searchsorted in numpy
                while (cog_res < NUM_SPECIES){
                    if (grain_mass_list(cog_res) >= (grain_mass_list(j) + grain_mass_list(k))){
                        break;
                    }
                    cog_res += 1;
                }
                if (j == k){
                    if (cog_res == NUM_SPECIES){
#pragma acc atomic update
                        Mmat[cog_res - 1] += 0.5 * (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res - 1) * numjtimesnumk;
                    }
                    else{
                        double eps = (grain_mass_list(j) + grain_mass_list(k) - grain_mass_list(cog_res - 1)) / (grain_mass_list(cog_res) - grain_mass_list(cog_res - 1));
#pragma acc atomic update
                        Mmat[cog_res]     += 0.5 * (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res) * numjtimesnumk * eps;
#pragma acc atomic update
                        Mmat[cog_res - 1] += 0.5 * (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res - 1) * numjtimesnumk * (1.0 - eps);
                    }
                }
                else{
                    if (cog_res == NUM_SPECIES){
#pragma acc atomic update
                        Mmat[cog_res - 1] += (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res - 1)* numjtimesnumk;
                    }
                    else{
                        double eps =  (grain_mass_list(j) + grain_mass_list(k) - grain_mass_list(cog_res - 1)) / (grain_mass_list(cog_res) - grain_mass_list(cog_res - 1));
#pragma acc atomic update
                        Mmat[cog_res]     += (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res)* numjtimesnumk * eps;
#pragma acc atomic update
                        Mmat[cog_res - 1] += (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res - 1)* numjtimesnumk * (1.0 - eps);
                    }
                }
                // gain via fragmentation
#pragma acc atomic update
                Mmat[0] += KL1[1] * grain_mass_list(k) / grain_mass_list(0) * numjtimesnumk;
#pragma acc atomic update
                Mmat[0] += KL2[1] * grain_mass_list(j) / grain_mass_list(0) * numjtimesnumk;
                // lost via coagulation and fragmentation
#pragma acc atomic update
                Mmat[k] -= (KL1[0] + KL1[1]) * numjtimesnumk;
                if (k != j){
#pragma acc atomic update
                    Mmat[j] -= (KL2[0] + KL2[1]) * numjtimesnumk;
                }
                // break if there is a problem
                /*
                if (isnan(Mmat(0))){
                }
                if (isnan(Mmat(j))){
                }
                if (isnan(Mmat(k))){
                }
                */
            }
        }
        dt_tot += dt_here;
#pragma acc loop vector
        for (int i = 0; i < NUM_SPECIES; ++ i){
            num_here[i] += Mmat[i] * dt_here;
            if (num_here[i] < 0)
            {
                num_here[i] = 0;
            }
        }
        //cout << dt_tot << '\t' << min_dt << '\t' << dt << '\n' << flush;
    }
#pragma acc loop vector
    for (int i = 0; i < NUM_SPECIES; i++){num[i] = num_here[i];}
}


void grain_growth(mesh &m, BootesArray<double> &stoppingtimemesh, double &dt){
    int NUMSPECIES = m.NUMSPECIES;
    //double *grain_number_array = new double[NUMSPECIES];
    //double *grain_vr_array     = new double[NUMSPECIES];
    //double *grain_vtheta_array = new double[NUMSPECIES];
    //double *grain_vphi_array   = new double[NUMSPECIES];
    double grain_number_array[NUMSPECIES];
    double grain_vr_array[NUMSPECIES];
    double grain_vtheta_array[NUMSPECIES];
    double grain_vphi_array[NUMSPECIES];
    double num_here[NUMSPECIES];
    double Mmat[NUMSPECIES];

	// #pragma omp parallel for collapse (3) schedule(dynamic) private (grain_number_array, grain_vr_array, grain_vtheta_array, grain_vphi_array)
	// add worker if need to pair vector length and work-load (in this case NUMSPECIES)
	//int NUM_GANGS = std::min((m.x3l-m.x3s)*(m.x2l-m.x2s)*(m.x1l-m.x1s), 1000);
    #pragma acc parallel loop gang worker collapse (3) vector_length(32) default (present) firstprivate(dt) \
        private (grain_number_array[0:NUMSPECIES], grain_vr_array[0:NUMSPECIES], grain_vtheta_array[0:NUMSPECIES], grain_vphi_array[0:NUMSPECIES], \
        num_here[0:NUMSPECIES], Mmat[0:NUMSPECIES])
	for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                /*
                if (rhomesh(j)[i] < 4e-17){         // if density is below some threshold, just skip coagulation calculations.
                    continue;
                }
                */

                #pragma acc loop vector
                for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
                    //double gas_rho = m.dcons(specIND, IDN, kk, jj, ii);
                    grain_number_array[specIND] = m.dcons(specIND, IDN, kk, jj, ii) / m.GrainMassList(specIND);
                    grain_vr_array[specIND]     = m.dcons(specIND, IM1, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                    grain_vtheta_array[specIND] = m.dcons(specIND, IM2, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                    grain_vphi_array[specIND]   = m.dcons(specIND, IM3, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                    //cout << m.dcons(specIND, IM1, kk, jj, ii) << '\t';
                }
                //cout << endl << flush;
                /** START OF GRAIN GROWTH IN ONE CELL **/
                /*
                double dt_here = dt;
            #pragma acc loop vector
                for (int i = 0; i < NUMSPECIES; i++){num_here[i] = grain_number_array[i];}
                bool redo = true;
                double dt_tot = 0;
                redo = false;   // set it to false first
            #pragma acc loop vector //collapse (2)
                for (int j = 0; j < NUMSPECIES; ++ j){
                    for (int k = j; k < NUMSPECIES; ++ k){
                        double dv_bulk = sqrt(pow((grain_vr_array[j] - grain_vr_array[k]), 2) + pow((grain_vtheta_array[j] - grain_vtheta_array[k]), 2) + pow((grain_vphi_array[j] - grain_vphi_array[k]), 2));
                        // double dv_vortex = dv_ormel(grain_size_list(j), grain_size_list(k), rhogas, tempgas);
                        double dv = dv_bulk;

                        double K = dv * M_PI * pow(m.GrainSizeList(j) + m.GrainSizeList(k), 2.0);

                        // replace the for statement here with new method, use V1.8's for if needed
                        // gain via coagulation
                        double numjtimesnumk = num_here[j] * num_here[k];
                        int cog_res = 0;        // the following loop is equivalent to searchsorted in numpy
                        while (cog_res < NUMSPECIES){
                            if (m.GrainMassList(cog_res) >= (m.GrainMassList(j) + m.GrainMassList(k))){
                                break;
                            }
                            cog_res += 1;
                        }
                        if (j == k){
                            if (cog_res == NUMSPECIES){
            #pragma acc atomic update
                                Mmat[cog_res - 1] += 0.5 * (K * m.GrainMassList(k) + K * m.GrainMassList(j)) / m.GrainMassList(cog_res - 1) * numjtimesnumk;
                            }
                            else{
                                double eps = (m.GrainMassList(j) + m.GrainMassList(k) - m.GrainMassList(cog_res - 1)) / (m.GrainMassList(cog_res) - m.GrainMassList(cog_res - 1));
            #pragma acc atomic update
                                Mmat[cog_res]     += 0.5 * (K * m.GrainMassList(k) + K * m.GrainMassList(j)) / m.GrainMassList(cog_res) * numjtimesnumk * eps;
            #pragma acc atomic update
                                Mmat[cog_res - 1] += 0.5 * (K * m.GrainMassList(k) + K * m.GrainMassList(j)) / m.GrainMassList(cog_res - 1) * numjtimesnumk * (1.0 - eps);
                            }
                        }
                        else{
                            if (cog_res == NUMSPECIES){
            #pragma acc atomic update
                                Mmat[cog_res - 1] += (K * m.GrainMassList(k) + K * m.GrainMassList(j)) / m.GrainMassList(cog_res - 1)* numjtimesnumk;
                            }
                            else{
                                double eps =  (m.GrainMassList(j) + m.GrainMassList(k) - m.GrainMassList(cog_res - 1)) / (m.GrainMassList(cog_res) - m.GrainMassList(cog_res - 1));
            #pragma acc atomic update
                                Mmat[cog_res]     += (K * m.GrainMassList(k) + K * m.GrainMassList(j)) / m.GrainMassList(cog_res)* numjtimesnumk * eps;
            #pragma acc atomic update
                                Mmat[cog_res - 1] += (K * m.GrainMassList(k) + K * m.GrainMassList(j)) / m.GrainMassList(cog_res - 1)* numjtimesnumk * (1.0 - eps);
                            }
                        }
                        // lost via coagulation
            #pragma acc atomic update
                        Mmat[k] -= K * numjtimesnumk;
                        if (k != j){
            #pragma acc atomic update
                            Mmat[j] -= K * numjtimesnumk;
                        }
                    }
                }
                dt_tot += dt_here;
            #pragma acc loop vector
                for (int i = 0; i < NUMSPECIES; ++ i){
                    num_here[i] += Mmat[i] * dt_here;
                    if (num_here[i] < 0)
                    {
                        num_here[i] = 0;
                    }
                }
                //cout << dt_tot << '\t' << min_dt << '\t' << dt << '\n' << flush;
            #pragma acc loop vector
                for (int i = 0; i < NUMSPECIES; i++){grain_number_array[i] = num_here[i];}
                */
                /** END OF GRAIN GROWTH IN ONE CELL **/
                grain_growth_one_cell_stick(grain_number_array,
                                      grain_vr_array, grain_vtheta_array, grain_vphi_array, num_here, Mmat,
                                      m.GrainSizeList, m.GrainMassList, dt, m.NUMSPECIES);

                // copy 1-cell results from grain_number_array to m.dcons

                #pragma acc loop vector
                for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++) {
                    if (grain_number_array[specIND] * m.GrainMassList(specIND) < m.dminDensity && false) {
                        m.dcons(specIND, IDN, kk, jj, ii) = m.dminDensity;
                        double rhogradphix1;
                        double rhogradphix2;
                        double rhogradphix3;
                        #ifdef ENABLE_GRAVITY
                        rhogradphix1 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x1surface(kk, jj, ii + 1) - m.grav->Phi_grav_x1surface(kk, jj, ii)) / m.dx1p(kk, jj, ii);
                        rhogradphix2 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x2surface(kk, jj + 1, ii) - m.grav->Phi_grav_x2surface(kk, jj, ii)) / m.dx2p(kk, jj, ii);
                        rhogradphix3 = m.dcons(specIND, IDN, kk, jj, ii) * (m.grav->Phi_grav_x3surface(kk + 1, jj, ii) - m.grav->Phi_grav_x3surface(kk, jj, ii)) / m.dx3p(kk, jj, ii);
                        #else   // set gravity to zero
                        rhogradphix1 = 0;
                        rhogradphix2 = 0;
                        rhogradphix3 = 0;
                        #endif // ENABLE_GRAVITY
                        #ifdef CARTESIAN_COORD
                        dust_terminalvelocityapprixmation_xyz(m.prim(IV1, kk, jj, ii), m.prim(IV2, kk, jj, ii), m.prim(IV3, kk, jj, ii),
                                                              rhogradphix1, rhogradphix2, rhogradphix3,
                                                              m.dcons(specIND, IDN, kk, jj, ii), stoppingtimemesh(specIND, kk, jj, ii),
                                                              m.dcons(specIND, IM1, kk, jj, ii), m.dcons(specIND, IM2, kk, jj, ii), m.dcons(specIND, IM3, kk, jj, ii)
                                                              );
                        #endif // CARTESIAN_COORD
                        #ifdef SPHERICAL_POLAR_COORD
                        dust_terminalvelocityapprixmation_rtp(m.prim(IV1, kk, jj, ii), m.prim(IV2, kk, jj, ii), m.prim(IV3, kk, jj, ii),
                                                              rhogradphix1, rhogradphix2, rhogradphix3,
                                                              m.dcons(specIND, IDN, kk, jj, ii), stoppingtimemesh(specIND, kk, jj, ii), m.x1v(ii), m.geo_cot(jj),
                                                              m.dcons(specIND, IM1, kk, jj, ii), m.dcons(specIND, IM2, kk, jj, ii), m.dcons(specIND, IM3, kk, jj, ii)
                                                              );
                        #endif // SPHERICAL_POLAR_COORD
                        #ifdef DEBUG
                        std::cout << "drho < 0:\t" << specIND << '\t' << kk << '\t' << jj << '\t' << ii << '\t'
                                  << m.dcons(specIND, IDN, kk, jj, ii) << '\t' << m.dcons(specIND, IM1, kk, jj, ii) << '\t'
                                  << m.dcons(specIND, IM2, kk, jj, ii) << '\t' << m.dcons(specIND, IM3, kk, jj, ii) << std::endl << flush;
                        #endif // DEBUG
                    }
                    else {
                        m.dcons(specIND, IDN, kk, jj, ii) = grain_number_array[specIND] * m.GrainMassList(specIND);
                        m.dcons(specIND, IM1, kk, jj, ii) = m.dcons(specIND, IDN, kk, jj, ii) * grain_vr_array[specIND];
                        m.dcons(specIND, IM2, kk, jj, ii) = m.dcons(specIND, IDN, kk, jj, ii) * grain_vtheta_array[specIND];
                        m.dcons(specIND, IM3, kk, jj, ii) = m.dcons(specIND, IDN, kk, jj, ii) * grain_vphi_array[specIND];
                    }
                }
            }
        }
    }
    /*
    delete[] grain_number_array;
    delete[] grain_vr_array;
    delete[] grain_vtheta_array;
    delete[] grain_vphi_array;
    */
    #pragma omp barrier
}
