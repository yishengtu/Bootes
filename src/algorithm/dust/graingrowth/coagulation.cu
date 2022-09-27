#include "coagulation.hpp"
#include "graingrowthmodel.hpp"
#include "../../util/util.hpp"
#include "../../index_def.hpp"
#include "../../mesh/mesh.hpp"
#include "../../../defs.hpp"
#include "../terminalvel.hpp"
#include <cmath>

__device__ void grain_growth_model_stick(double &s1, double &s2, double &dv, double res[2]){
    res[0] = dv * M_PI * pow((s1 + s2), 2.0);
    res[1] = 0.0;
}


__global__ void grain_growth_one_cell(double *num,
                           double *vr, double *vtheta, double *vphi, double *num_here, double *Mmat,
                           double *grain_size_list, double *grain_mass_list, double dt, int NUM_SPECIES){
    double dt_here = dt;

    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < NUM_SPECIES; i+=stride) {num_here[i] = num[i];}
    __syncthreads();

    bool redo = true;
    double dt_tot = 0;
    while (redo || dt_tot < dt)
    {
        redo = false;   // set it to false first
        // BootesArray<double> Mmat; Mmat.NewBootesArray(NUM_SPECIES);         // initialize with 0s

        for (int k = index; k < NUM_SPECIES; k+=stride){
            Mmat[k] = 0;
        }
	  __syncthreads();

        for (int j = index; j < NUM_SPECIES; j+=stride){
            for (int k = j; k < NUM_SPECIES; ++ k){
		//if (k < j) continue;
                double dv_bulk = sqrt(pow((vr[j] - vr[k]), 2) + pow((vtheta[j] - vtheta[k]), 2) + pow((vphi[j] - vphi[k]), 2));
                // double dv_vortex = dv_ormel(grain_size_list[j], grain_size_list[k], rhogas, tempgas);
                double dv = dv_bulk;

                double KL1[2];                  //K1 = KL1[0], L1 = KL1[1]
                double KL2[2];
                grain_growth_model_stick(grain_size_list[j], grain_size_list[k], dv, KL1);
                grain_growth_model_stick(grain_size_list[k], grain_size_list[j], dv, KL2);

                // replace the for statement here with new method, use V1.8's for if needed
                // gain via coagulation
                double numjtimesnumk = num_here[j] * num_here[k];
                int cog_res = 0;        // the following loop is equivalent to searchsorted in numpy
                while (cog_res < NUM_SPECIES){
                    if (grain_mass_list[cog_res] >= (grain_mass_list[j] + grain_mass_list[k])){
                        break;
                    }
                    cog_res += 1;
                }
                if (j == k){
                    if (cog_res == NUM_SPECIES){
                        Mmat[cog_res - 1] += 0.5 * (KL1[0] * grain_mass_list[k] + KL2[0] * grain_mass_list[j]) / grain_mass_list[cog_res - 1] * numjtimesnumk;
                    }
                    else{
                        double eps = (grain_mass_list[j] + grain_mass_list[k] - grain_mass_list[cog_res - 1]) / (grain_mass_list[cog_res] - grain_mass_list[cog_res - 1]);
                        Mmat[cog_res]     += 0.5 * (KL1[0] * grain_mass_list[k] + KL2[0] * grain_mass_list[j]) / grain_mass_list[cog_res] * numjtimesnumk * eps;
                        Mmat[cog_res - 1] += 0.5 * (KL1[0] * grain_mass_list[k] + KL2[0] * grain_mass_list[j]) / grain_mass_list[cog_res - 1] * numjtimesnumk * (1.0 - eps);
                    }
                }
                else{
                    if (cog_res == NUM_SPECIES){
                        Mmat[cog_res - 1] += (KL1[0] * grain_mass_list[k] + KL2[0] * grain_mass_list[j]) / grain_mass_list[cog_res - 1]* numjtimesnumk;
                    }
                    else{
                        double eps =  (grain_mass_list[j] + grain_mass_list[k] - grain_mass_list[cog_res - 1]) / (grain_mass_list[cog_res] - grain_mass_list[cog_res - 1]);
                        Mmat[cog_res]     += (KL1[0] * grain_mass_list[k] + KL2[0] * grain_mass_list[j]) / grain_mass_list[cog_res]* numjtimesnumk * eps;
                        Mmat[cog_res - 1] += (KL1[0] * grain_mass_list[k] + KL2[0] * grain_mass_list[j]) / grain_mass_list[cog_res - 1]* numjtimesnumk * (1.0 - eps);
                    }
                }
                // gain via fragmentation
                Mmat[0] += KL1[1] * grain_mass_list[k] / grain_mass_list[0] * numjtimesnumk;
                Mmat[0] += KL2[1] * grain_mass_list[j] / grain_mass_list[0] * numjtimesnumk;
                // lost via coagulation and fragmentation
                Mmat[k] -= (KL1[0] + KL1[1]) * numjtimesnumk;
                if (k != j){
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
        __syncthreads();

        dt_tot += dt_here;
        for (int i = index; i < NUM_SPECIES; i+=stride){
            num_here[i] += Mmat[i] * dt_here;
            if (num_here[i] < 0)
            {
                num_here[i] = 0;
            }
        }
	  __syncthreads();

        //cout << dt_tot << '\t' << min_dt << '\t' << dt << '\n' << flush;
    }
      __syncthreads();

    for (int i = index; i < NUM_SPECIES; i+=stride){num[i] = num_here[i];}
}


void grain_growth(mesh &m, BootesArray<double> &stoppingtimemesh, double &dt){
    int NUMSPECIES = m.NUMSPECIES;
    //double *grain_number_array = new double[NUMSPECIES];
    //double *grain_vr_array     = new double[NUMSPECIES];
    //double *grain_vtheta_array = new double[NUMSPECIES];
    //double *grain_vphi_array   = new double[NUMSPECIES];
    double *grain_number_array;
    double *grain_vr_array;
    double *grain_vtheta_array;
    double *grain_vphi_array;
    double *num_here;
    double *Mmat;
    

    size_t nspbytes = NUMSPECIES*sizeof(double);
    
    int device_id;
    
    cudaGetDevice(&device_id);

    cudaMallocManaged(&grain_number_array, nspbytes);
    cudaMallocManaged(&grain_vr_array, nspbytes);
    cudaMallocManaged(&grain_vtheta_array, nspbytes);
    cudaMallocManaged(&grain_vphi_array, nspbytes);
    cudaMallocManaged(&num_here, nspbytes);
    cudaMallocManaged(&Mmat, nspbytes);

	for (int kk = m.x3s; kk < m.x3l; kk ++){
        for (int jj = m.x2s; jj < m.x2l; jj ++){
            for (int ii = m.x1s; ii < m.x1l; ii ++){
                /*
                if (rhomesh(j)[i] < 4e-17){         // if density is below some threshold, just skip coagulation calculations.
                    continue;
                }
                */
                for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
                    //double gas_rho = m.dcons(specIND, IDN, kk, jj, ii);
                    grain_number_array[specIND] = m.dcons(specIND, IDN, kk, jj, ii) / m.GrainMassList(specIND);
                    grain_vr_array[specIND]     = m.dcons(specIND, IM1, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                    grain_vtheta_array[specIND] = m.dcons(specIND, IM2, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                    grain_vphi_array[specIND]   = m.dcons(specIND, IM3, kk, jj, ii) / m.dcons(specIND, IDN, kk, jj, ii);
                    //cout << m.dcons(specIND, IM1, kk, jj, ii) << '\t';
                }
		    cudaMemPrefetchAsync(grain_number_array, nspbytes, device_id);
		    cudaMemPrefetchAsync(grain_vr_array, nspbytes, device_id);
		    cudaMemPrefetchAsync(grain_vtheta_array, nspbytes, device_id);
		    cudaMemPrefetchAsync(grain_vphi_array, nspbytes, device_id);
		    cudaMemPrefetchAsync(num_here, nspbytes, device_id);
		    cudaMemPrefetchAsync(Mmat, nspbytes, device_id);
                //cout << endl << flush;
                size_t grnbytes = m.GrainSizeList.shape()[0]*sizeof(double);
		double *d_GrainSizeList_arr;
		double *d_GrainMassList_arr;
		cudaMalloc(&d_GrainSizeList_arr, grnbytes);
		cudaMalloc(&d_GrainMassList_arr, grnbytes);
                cudaMemcpy(d_GrainSizeList_arr, m.GrainSizeList.data(), grnbytes, cudaMemcpyHostToDevice);
                cudaMemcpy(d_GrainMassList_arr, m.GrainMassList.data(), grnbytes, cudaMemcpyHostToDevice);
                int MB=1;
                grain_growth_one_cell<<<NUMSPECIES/MB, MB>>>(grain_number_array,
                                      grain_vr_array, grain_vtheta_array, grain_vphi_array, num_here, Mmat,
                                      d_GrainSizeList_arr, d_GrainMassList_arr, dt, m.NUMSPECIES);
                // copy 1-cell results from grain_number_array to m.dcons

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
    /*
    delete[] grain_number_array;
    delete[] grain_vr_array;
    delete[] grain_vtheta_array;
    delete[] grain_vphi_array;
    */
}
