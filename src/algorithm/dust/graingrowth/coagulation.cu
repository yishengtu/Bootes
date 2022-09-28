#include "coagulation.hpp"
#include "graingrowthmodel.hpp"
#include "../../util/util.hpp"
#include "../../index_def.hpp"
#include "../../mesh/mesh.hpp"
#include "../../../defs.hpp"
//#include "../terminalvel.hpp"
#include <cmath>

#define CHECK_CUDA_ERROR(val) check( (val), #val, __FILE__, __LINE__)
template <class T>
void check(T err, const char * const errStr, const char * const file, const int line) {
    if(err != cudaSuccess) {
	    cerr << "CUDA error at: " << file << ":" << line << endl;
	    cerr << cudaGetErrorString(err) << " : " << errStr << endl;
	    exit(1);
    }
}

struct meshsim
{
double *Phi_grav_x1surface;
double *Phi_grav_x2surface;
double *Phi_grav_x3surface;

double *dx1p;
double *dx2p;
double *dx3p;
meshsim(int size)
{
   int state = 0;
   cudaMalloc(&Phi_grav_x1surface, size*sizeof(double));
   cudaMalloc(&Phi_grav_x2surface, size*sizeof(double));
   cudaMalloc(&Phi_grav_x3surface, size*sizeof(double));
   cudaMalloc(&dx1p, size*sizeof(double));
   cudaMalloc(&dx2p, size*sizeof(double));
   cudaMalloc(&dx3p, size*sizeof(double));
   //cudaMemcpy(Phi_grav_x1surface, &state, sizeof(double), cudaMemcpyHostToDevice);
}

~meshsim()
{
   cudaFree(Phi_grav_x1surface);
   cudaFree(Phi_grav_x2surface);
   cudaFree(Phi_grav_x3surface);
   cudaFree(dx1p);
   cudaFree(dx2p);
   cudaFree(dx3p);
}

/*
__device__ void lock() {
   while(atomicCAS(mutex,0,1) != 0);
}

__device__ void unlock() {
   atomicExch(mutex, 0);
}
*/
};


__device__ void grain_growth_model_stick(double &s1, double &s2, double &dv, double res[2]){
    res[0] = dv * M_PI * pow((s1 + s2), 2.0);
    res[1] = 0.0;
}

__device__ void cu_dust_terminalvelocityapprixmation_xyz(double &vg1, double &vg2, double &vg3,
                                           double &g1,  double &g2,  double &g3,
                                           double &rhod, double &ts,
                                           double &pd1, double &pd2, double &pd3){
    pd1 = rhod * vg1 + g1 * ts;
    pd2 = rhod * vg2 + g2 * ts;
    pd3 = rhod * vg3 + g3 * ts;
    }

__global__ void growth(double *dcons, double *prim, double *stoppingtimemesh, meshsim grav, double dt, double NUMSPECIES, double dminDensity, int *shape,
		int x1s, int x1l, int x2s, int x2l, int x3s, int x3l,
                //double *grain_number_array, double *grain_vr_array, double *grain_vtheta_array,
                //double *grain_vphi_array, double *num_here, double *Mmat,
		double *d_GrainSizeList_arr, double *d_GrainMassList_arr);

__device__ void grain_growth_one_cell(double *num,
                           double *vr, double *vtheta, double *vphi, double *num_here, double *Mmat,
                           double *grain_size_list, double *grain_mass_list, double dt, int NUM_SPECIES){
    double dt_here = dt;
    double temp;
    int index = threadIdx.x;// + blockIdx.x * blockDim.x;
    int stride = blockDim.x;// * gridDim.x;

    for (int i = index; i < NUM_SPECIES; i+=stride) {num_here[i] = num[i];}
    //__syncthreads();

    bool redo = true;
    double dt_tot = 0;
    while (redo || dt_tot < dt)
    {
        redo = false;   // set it to false first
        // BootesArray<double> Mmat; Mmat.NewBootesArray(NUM_SPECIES);         // initialize with 0s

        for (int k = index; k < NUM_SPECIES; k+=stride){
            Mmat[k] = 0;
        }
	  //__syncthreads();

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
                        temp = (KL1[0] * grain_mass_list[k] + KL2[0] * grain_mass_list[j]) / grain_mass_list[cog_res - 1]* numjtimesnumk;
			//Mmat[cog_res - 1] += temp;
			atomicAdd ( Mmat[cog_res - 1],temp);
                    //    Mmat[cog_res - 1] += (KL1[0] * grain_mass_list[k] + KL2[0] * grain_mass_list[j]) / grain_mass_list[cog_res - 1]* numjtimesnumk;
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
        //__syncthreads();

        dt_tot += dt_here;
        for (int i = index; i < NUM_SPECIES; i+=stride){
            num_here[i] += Mmat[i] * dt_here;
            if (num_here[i] < 0)
            {
                num_here[i] = 0;
            }
        }
	  //__syncthreads();

        //cout << dt_tot << '\t' << min_dt << '\t' << dt << '\n' << flush;
    }
      //__syncthreads();

    for (int i = index; i < NUM_SPECIES; i+=stride){num[i] = num_here[i];}
}


void grain_growth(mesh &m, BootesArray<double> &stoppingtimemesh, double &dt){
    int NUMSPECIES = m.NUMSPECIES;
    double dminDensity = m.dminDensity;
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
    //cudaError_t cerr;

    size_t nspbytes = NUMSPECIES*sizeof(double);
    int device_id;
    
    CHECK_CUDA_ERROR(cudaGetDevice(&device_id));

    CHECK_CUDA_ERROR(cudaMallocManaged(&grain_number_array, nspbytes));
    //std::cout<<"cerr "<<cerr<<std::endl;
    CHECK_CUDA_ERROR(cudaMallocManaged(&grain_vr_array, nspbytes));
    CHECK_CUDA_ERROR(cudaMallocManaged(&grain_vtheta_array, nspbytes));
    CHECK_CUDA_ERROR(cudaMallocManaged(&grain_vphi_array, nspbytes));
    CHECK_CUDA_ERROR(cudaMallocManaged(&num_here, nspbytes));
    CHECK_CUDA_ERROR(cudaMallocManaged(&Mmat, nspbytes));
    
    int NG=2; 
    int size1=m.dcons.shape()[4]-2*NG, size2=m.dcons.shape()[3]-2*NG, size3=m.dcons.shape()[2]-2*NG;
    int ncell = m.dcons.shape()[4]*m.dcons.shape()[3]*m.dcons.shape()[2];

    //std::cout<<"size1 "<<size1<<" size2 "<<size2<<" size3 "<<size3<<std::endl;
    int BLKX=32, BLKY=32, BLKZ=32;
    int MB=32;
    size_t grnbytes = m.GrainSizeList.shape()[0]*sizeof(double);
    double *d_GrainSizeList_arr;
    double *d_GrainMassList_arr;
    double *d_stoppingtimemesh;
    size_t stpbytes = stoppingtimemesh.shape()[0]*stoppingtimemesh.shape()[1]*stoppingtimemesh.shape()[2]*stoppingtimemesh.shape()[3]*sizeof(double);
    CHECK_CUDA_ERROR(cudaMalloc(&d_GrainSizeList_arr, grnbytes));
    CHECK_CUDA_ERROR(cudaMalloc(&d_GrainMassList_arr, grnbytes));
    CHECK_CUDA_ERROR(cudaMalloc(&d_stoppingtimemesh, stpbytes));
    CHECK_CUDA_ERROR(cudaMemcpy(d_GrainSizeList_arr, m.GrainSizeList.data(), grnbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(d_GrainMassList_arr, m.GrainMassList.data(), grnbytes, cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(d_stoppingtimemesh, stoppingtimemesh.data(), stpbytes, cudaMemcpyHostToDevice));

    meshsim grav=meshsim(ncell);
    CHECK_CUDA_ERROR(cudaMemcpy(grav.Phi_grav_x1surface,m.grav->Phi_grav_x1surface.data(),ncell*sizeof(double),cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(grav.Phi_grav_x2surface,m.grav->Phi_grav_x2surface.data(),ncell*sizeof(double),cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(grav.Phi_grav_x3surface,m.grav->Phi_grav_x3surface.data(),ncell*sizeof(double),cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(grav.dx1p,m.dx1p.data(),ncell*sizeof(double),cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(grav.dx2p,m.dx2p.data(),ncell*sizeof(double),cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(grav.dx3p,m.dx3p.data(),ncell*sizeof(double),cudaMemcpyHostToDevice));
    

    double *d_dcons = (double*)acc_deviceptr(m.dcons.data());
    double *d_prim = (double*)acc_deviceptr(m.prim.data());

    int *shape_m = m.dcons.shape();
    dim3 BlocksperGrid(size1/BLKX+1,size2/BLKY+1, size3/BLKZ+1);
    dim3 ThreadsperBlock(BLKX,BLKY,BLKZ); 
    int x1s = m.x1s;int x1l = m.x1l;int x2s = m.x2s;int x2l = m.x2l;int x3s = m.x3s;int x3l = m.x3l;
    int sharedbytes = 6*NUMSPECIES*sizeof(double); 
    //std::cout<<"x1s "<<x1s<<" x1l "<<x1l<<" x2s "<<x2s<<" x2l "<<x2l<<" x3s "<<x3s<<" x3l "<<x3l<<std::endl;
    growth<<<BlocksperGrid, ThreadsperBlock, sharedbytes>>>(d_dcons, d_prim, d_stoppingtimemesh, grav, dt, NUMSPECIES, dminDensity, shape_m, x1s, x1l, x2s, x2l, x3s, x3l, d_GrainSizeList_arr, d_GrainMassList_arr);
    cudaDeviceSynchronize();

}

__global__ void growth(double *dcons, double *prim, double *stoppingtimemesh, meshsim grav, double dt, double NUMSPECIES,double dminDensity, int *shape,
		                int x1s, int x1l, int x2s, int x2l, int x3s, int x3l,
		//double *grain_number_array, double *grain_vr_array, double *grain_vtheta_array, 
		//double *grain_vphi_array, double *num_here, double *Mmat,
		double *GrainSizeList_arr, double *GrainMassList_arr){

    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    //int index = blockIdx.x;
    int stridex = blockDim.x * gridDim.x;
    int stridey = blockDim.y * gridDim.y;
    int stridez = blockDim.z * gridDim.z;
    int index_x = threadIdx.x + blockDim.x*blockIdx.x;
    int index_y = threadIdx.y + blockDim.y*blockIdx.y;
    int index_z = threadIdx.z + blockDim.z*blockIdx.z;
    //size_t nspbytes = NUMSPECIES*sizeof(double);
    
    index_x = blockId.x;
    index_y = blockId.y;
    index_z = blockId.z;

    stridex = gridDim.x;
    stridey = gridDim.y;
    stridez = gridDim.z;

    int device_id;
    int size1 = shape[4];
    int size2 = shape[3];
    int size3 = shape[2];
    int size4 = shape[1];
    
    __shared__ grain_number_array[];
    __shared__ grain_vr_array[];
    __shared__ grain_vtheta_array[];
    __shared__ grain_vphi_array[];
    __shared__ num_here[];
    __shared__ Mmat[];

    //cudaGetDevice(&device_id);

    for (int kk = index_z+x3s; kk < x3l; kk +=stridez){
        for (int jj = index_y+x2s; jj <x2l; jj +=stridey){
            for (int ii = index_x+x1s; ii < x1l; ii +=stridex){
    
    /*int kk = index_z+x3s;
    int jj = index_y+x2s;
    int ii = index_x+x1s;
    if (kk >= x3l) return;
    if (jj >= x2l) return;
    if (ii >= x1l) return;
    */
                /*
                if (rhomesh(j)[i] < 4e-17){         // if density is below some threshold, just skip coagulation calculations.
                    continue;
                }
                */
                for (int specIND = theadId.x; specIND < NUMSPECIES; specIND +=blockDim.x){
                    //double gas_rho = m.dcons(specIND, IDN, kk, jj, ii);
	            //std::cout<<"grain_number_array "<<grain_number_array[specIND]<<std::endl;
		    int idx_IDN = ii + size1*(jj + size2*(kk + size3 * (IDN + size4 * specIND)));
		    int idx_IM1 = ii + size1*(jj + size2*(kk + size3 * (IM1 + size4 * specIND)));
		    int idx_IM2 = ii + size1*(jj + size2*(kk + size3 * (IM2 + size4 * specIND)));
		    int idx_IM3 = ii + size1*(jj + size2*(kk + size3 * (IM3 + size4 * specIND)));
                    grain_number_array[specIND] = dcons[idx_IDN] / GrainMassList_arr[specIND];
                    grain_vr_array[specIND]     = dcons[idx_IM1] / dcons[idx_IDN];
                    grain_vtheta_array[specIND] = dcons[idx_IM2] / dcons[idx_IDN];
                    grain_vphi_array[specIND]   = dcons[idx_IM3] / dcons[idx_IDN];
		    
                    //cout << m.dcons(specIND, IM1, kk, jj, ii) << '\t';
                }
                __syncthreads();
		//    cudaMemPrefetchAsync(grain_number_array, nspbytes, device_id);
		//    cudaMemPrefetchAsync(grain_vr_array, nspbytes, device_id);
		//    cudaMemPrefetchAsync(grain_vtheta_array, nspbytes, device_id);
		//    cudaMemPrefetchAsync(grain_vphi_array, nspbytes, device_id);
		//    cudaMemPrefetchAsync(num_here, nspbytes, device_id);
		//    cudaMemPrefetchAsync(Mmat, nspbytes, device_id);
                //cout << endl << flush;

              
                grain_growth_one_cell(grain_number_array,
                                      grain_vr_array, grain_vtheta_array, grain_vphi_array, num_here, Mmat,
                                      GrainSizeList_arr, GrainMassList_arr, dt, NUMSPECIES);
		__syncthreads();
                // copy 1-cell results from grain_number_array to m.dcons

                for (int specIND = 0; specIND < NUMSPECIES; specIND ++) {
		    int idx_IDN = ii + size1*(jj + size2*(kk + size3 * (IDN + size4 * specIND)));
		    int idx_IM1 = ii + size1*(jj + size2*(kk + size3 * (IM1 + size4 * specIND)));
		    int idx_IM2 = ii + size1*(jj + size2*(kk + size3 * (IM2 + size4 * specIND)));
		    int idx_IM3 = ii + size1*(jj + size2*(kk + size3 * (IM3 + size4 * specIND)));

                    int idx_IV1 = ii + size1*(jj + size2*(kk + size3 * IV1)); 
                    int idx_IV2 = ii + size1*(jj + size2*(kk + size3 * IV2)); 
                    int idx_IV3 = ii + size1*(jj + size2*(kk + size3 * IV3)); 
		    
		    int idx_stop= ii + size1*(jj + size2*(kk + size3 * specIND));

		    if (grain_number_array[specIND] * GrainMassList_arr[specIND] < dminDensity) {
                        dcons[idx_IDN] = dminDensity;
                        double rhogradphix1;
                        double rhogradphix2;
                        double rhogradphix3;
                        #ifdef ENABLE_GRAVITY
			int kji = ii + size1*(jj + size2*kk);
			int kji1 = ii + 1 + size1*(jj + size2*kk);
			int kj1i = ii + size1*(jj + 1 + size2*kk);
			int k1ji = ii + size1*(jj + (size2*kk+1));

                        rhogradphix1 = dcons[idx_IDN] * (grav.Phi_grav_x1surface[kji1] - grav.Phi_grav_x1surface[kji]) / grav.dx1p[kji];
                        rhogradphix2 = dcons[idx_IDN] * (grav.Phi_grav_x2surface[kj1i] - grav.Phi_grav_x2surface[kji]) / grav.dx2p[kji];
                        rhogradphix3 = dcons[idx_IDN] * (grav.Phi_grav_x3surface[k1ji] - grav.Phi_grav_x3surface[kji]) / grav.dx3p[kji];
                        //rhogradphix1 = dcons[idx_IDN] * (m.grav->Phi_grav_x1surface(kk, jj, ii + 1) - m.grav->Phi_grav_x1surface(kk, jj, ii)) / m.dx1p(kk, jj, ii);
                        //rhogradphix2 = dcons[idx_IDN] * (m.grav->Phi_grav_x2surface(kk, jj + 1, ii) - m.grav->Phi_grav_x2surface(kk, jj, ii)) / m.dx2p(kk, jj, ii);
                        //rhogradphix3 = dcons[idx_IDN] * (m.grav->Phi_grav_x3surface(kk + 1, jj, ii) - m.grav->Phi_grav_x3surface(kk, jj, ii)) / m.dx3p(kk, jj, ii);
                        #else   // set gravity to zero
                        rhogradphix1 = 0;
                        rhogradphix2 = 0;
                        rhogradphix3 = 0;
                        #endif // ENABLE_GRAVITY
                        #ifdef CARTESIAN_COORD
                        cu_dust_terminalvelocityapprixmation_xyz(prim[idx_IV1], prim[idx_IV2], prim[idx_IV3],
                                                              rhogradphix1, rhogradphix2, rhogradphix3,
                                                              dcons[idx_IDN], stoppingtimemesh[idx_stop],
                                                              dcons[idx_IM1], dcons[idx_IM2], dcons[idx_IM3]
                                                              );
                        #endif // CARTESIAN_COORD
                        #ifdef SPHERICAL_POLAR_COORD
                        dust_terminalvelocityapprixmation_rtp(m.prim(IV1, kk, jj, ii), m.prim(IV2, kk, jj, ii), m.prim(IV3, kk, jj, ii),
                                                              rhogradphix1, rhogradphix2, rhogradphix3,
                                                              m.dcons(specIND, IDN, kk, jj, ii), stoppingtimemesh(specIND, kk, jj, ii), m.x1v(ii), m.geo_cot(jj),
                                                              m.dcons(specIND, IM1, kk, jj, ii), m.dcons(specIND, IM2, kk, jj, ii), m.dcons(specIND, IM3, kk, jj, ii)
                                                              );
                        #endif // SPHERICAL_POLAR_COORD
                        //#ifdef DEBUG
                        //std::cout << "drho < 0:\t" << specIND << '\t' << kk << '\t' << jj << '\t' << ii << '\t'
                        //          << m.dcons(specIND, IDN, kk, jj, ii) << '\t' << m.dcons(specIND, IM1, kk, jj, ii) << '\t'
                        //          << m.dcons(specIND, IM2, kk, jj, ii) << '\t' << m.dcons(specIND, IM3, kk, jj, ii) << std::endl << flush;
                        //#endif // DEBUG
                    }
                    else {
                        dcons[idx_IDN] = grain_number_array[specIND] * GrainMassList_arr[specIND];
                        dcons[idx_IM1] = dcons[idx_IDN] * grain_vr_array[specIND];
                        dcons[idx_IM2] = dcons[idx_IDN] * grain_vtheta_array[specIND];
                        dcons[idx_IM3] = dcons[idx_IDN] * grain_vphi_array[specIND];
                    }
                }	
    //        }
    //    }
   // }
    /*
    delete[] grain_number_array;
    delete[] grain_vr_array;
    delete[] grain_vtheta_array;
    delete[] grain_vphi_array;
    */
}
