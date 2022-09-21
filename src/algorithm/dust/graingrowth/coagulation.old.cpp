#ifdef SKIP
#include "coagulation.hpp"
#include "graingrowthmodel.hpp"
#include "../../util/util.hpp"
#include "../../index_def.hpp"
#include "../../mesh/mesh.hpp"
#include "../../../defs.hpp"
#include "../terminalvel.hpp"
#include <cmath>


#pragma acc routine seq
void grain_growth_one_cell(double* num,
                           double* vr, double* vtheta, double* vphi,
                           BootesArray<double> &grain_size_list, BootesArray<double> &grain_mass_list,
                           int NUM_SPECIES, double dt){
    // TODO: parallel this function
    double dt_here = dt;
    double num_here[NUM_SPECIES];
    for (int i = 0; i < NUM_SPECIES; i++){num_here[i] = num[i];}
    bool redo = true;

    double Mmat[NUM_SPECIES];         // initialize with 0s
    for (int k = 0; k < NUM_SPECIES; k++){
        Mmat[k] = 0;
    }
    for (int j = 0; j < NUM_SPECIES; j++){
        for (int k = j; k < NUM_SPECIES; k++){
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
            // TODO: cog_res should start from 0, but 0 giving a error.
            int cog_res = 0;        // the following loop is equivalent to searchsorted in numpy
            while (cog_res < NUM_SPECIES){
                if (grain_mass_list(cog_res) >= (grain_mass_list(j) + grain_mass_list(k))){
                    break;
                }
                cog_res += 1;
            }
            if (j == k){
                if (cog_res == NUM_SPECIES){
                    Mmat[cog_res - 1] += 0.5 * (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res - 1) * numjtimesnumk;
                }
                else{
                    double eps = (grain_mass_list(j) + grain_mass_list(k) - grain_mass_list(cog_res - 1)) / (grain_mass_list(cog_res) - grain_mass_list(cog_res - 1));
                    Mmat[cog_res]     += 0.5 * (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res) * numjtimesnumk * eps;
                    Mmat[cog_res - 1] += 0.5 * (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res - 1) * numjtimesnumk * (1.0 - eps);
                }
            }
            else{
                if (cog_res == NUM_SPECIES){
                    Mmat[cog_res - 1] += (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res - 1) * numjtimesnumk;
                }
                else{
                    double eps =  (grain_mass_list(j) + grain_mass_list(k) - grain_mass_list(cog_res - 1)) / (grain_mass_list(cog_res) - grain_mass_list(cog_res - 1));
                    Mmat[cog_res]     += (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res)* numjtimesnumk * eps;
                    Mmat[cog_res - 1] += (KL1[0] * grain_mass_list(k) + KL2[0] * grain_mass_list(j)) / grain_mass_list(cog_res - 1)* numjtimesnumk * (1.0 - eps);
                }
            }
            // gain via fragmentation
            Mmat[0] += KL1[1] * grain_mass_list(k) / grain_mass_list(0) * numjtimesnumk;
            Mmat[0] += KL2[1] * grain_mass_list(j) / grain_mass_list(0) * numjtimesnumk;
            // lost via coagulation and fragmentation
            Mmat[k] -= (KL1[0] + KL1[1]) * numjtimesnumk;
            if (k != j){
                Mmat[j] -= (KL2[0] + KL2[1]) * numjtimesnumk;
            }
            #ifndef _OPENACC
            // replacement complete
            if (isnan(Mmat[0])){
                cout << "0 problem \n";
                cout << KL1[0] << '\t' << KL1[1] << '\t' << KL2[0] << '\t' << KL2[1] << '\t' << grain_size_list(j) << '\t' << grain_size_list(k) << '\t'
                     << vr[j] << '\t' << vr[k] << '\t' <<  vtheta[j] << '\t' << vtheta[k] << '\t' <<  vphi[j] << '\t' << vphi[k] << '\n';
                throw 1;
            }
            if (isnan(Mmat[j])){
                cout << "j problem \n";
                cout << KL1[0] << '\t' << KL1[1] << '\t' << KL2[0] << '\t' << KL2[1] << '\t' << grain_size_list(j) << '\t' << grain_size_list(k) << '\t' << dv << '\n';
                throw 1;
            }
            if (isnan(Mmat[k])){
                cout << "k problem \n";
                cout << KL1[0] << '\t' << KL1[1] << '\t' << KL2[0] << '\t' << KL2[1] << '\t' << grain_size_list(j) << '\t' << grain_size_list(k) << '\t' << dv << '\n';
                throw 1;
            }
            #endif // _OPENACC
        }
    }
    for (int i = 0; i < NUM_SPECIES; i ++){
        num_here[i] += Mmat[i] * dt_here;
        if (num_here[i] < 0)
        {
            num_here[i] = 0;
        }
    }
    #ifndef _OPENACC
    for (int i = 0; i < NUM_SPECIES; ++ i){
        if (isnan(num_here[i])){
            for (int i = 0; i < NUM_SPECIES; ++ i){
                cout << num[i] << '\t';
            }
            cout << '\n';
            for (int i = 0; i < NUM_SPECIES; ++ i){
                cout << Mmat[i] << '\t';
            }
            cout << '\n';
            cout << "END END END" << '\n';
            throw 1;
        }
    }
    #endif // _OPENACC
    for (int i = 0; i < NUM_SPECIES; i++){num[i] = num_here[i];}
}


void grain_growth(mesh &m, BootesArray<double> &stoppingtimemesh, double &dt){
	//#pragma omp parallel for collapse (3) schedule(dynamic)
    /*
    double grain_number_array[m.NUMSPECIES];
    double grain_vr_array[m.NUMSPECIES];
    double grain_vtheta_array[m.NUMSPECIES];
    double grain_vphi_array[m.NUMSPECIES];
    double num_here[m.NUMSPECIES];
    double Mmat[m.NUMSPECIES];         // initialize with 0s
    double GrainSizeList[m.NUMSPECIES];
    double GrainMassList[m.NUMSPECIES];
    for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
        GrainMassList[specIND] = m.GrainMassList(specIND);
        GrainSizeList[specIND] = m.GrainSizeList(specIND);
    }
    */
    int NUMSPECIES = m.NUMSPECIES;
    /*#pragma acc parallel loop gang collapse (3) default (present) private (grain_number_array[0:NUMSPECIES], \
         grain_vr_array[0:NUMSPECIES], grain_vtheta_array[0:NUMSPECIES], grain_vphi_array[0:NUMSPECIES], num_here[0:NUMSPECIES], Mmat[0:NUMSPECIES], GrainSizeList[0:NUMSPECIES], GrainMassList[0:NUMSPECIES]) \
         firstprivate(dt)
    */
    #pragma acc parallel loop gang collapse (3) default (present) firstprivate(dt)
	for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                /*
                if (rhomesh[j][i] < 4e-17){         // if density is below some threshold, just skip coagulation calculations.
                    continue;
                }
                */
                double grain_number_array[m.NUMSPECIES];
                double grain_vr_array[m.NUMSPECIES];
                double grain_vtheta_array[m.NUMSPECIES];
                double grain_vphi_array[m.NUMSPECIES];
                #pragma acc loop vector
                for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
                    double gas_rho = m.dcons(specIND, IDN, kk, jj, ii);
                    grain_number_array[specIND] = gas_rho / m.GrainMassList(specIND);
                    grain_vr_array[specIND]     = m.dcons(specIND, IM1, kk, jj, ii) / gas_rho;
                    grain_vtheta_array[specIND] = m.dcons(specIND, IM2, kk, jj, ii) / gas_rho;
                    grain_vphi_array[specIND]   = m.dcons(specIND, IM3, kk, jj, ii) / gas_rho;
                    //cout << m.dcons(specIND, IM1, kk, jj, ii) << '\t';
                }

                //cout << endl << flush;
                grain_growth_one_cell(grain_number_array,
                                     grain_vr_array, grain_vtheta_array, grain_vphi_array,
                                     m.GrainSizeList, m.GrainMassList, m.NUMSPECIES, dt);
                // copy 1-cell results from grain_number_array to m.dcons

                #pragma acc loop vector
                for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++) {
                    if (grain_number_array[specIND] * m.GrainMassList(specIND) < m.dminDensity) {
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
    #pragma omp barrier
}
#endif
