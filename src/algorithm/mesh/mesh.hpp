#ifndef MESH_HPP_
#define MESH_HPP_
#include "../BootesArray.hpp"
#include "../gravity/gravity.hpp"
#include "../physical_constants.hpp"


class mesh{
    public:
        /** consts **/
        PhysicalConst pconst;
        /** grid **/
        int dim;
        BootesArray<double> x1v;       // cell center (1D array)
        BootesArray<double> x2v;
        BootesArray<double> x3v;
        BootesArray<double> x1f;       // cell face (1D array)
        BootesArray<double> x2f;
        BootesArray<double> x3f;
        BootesArray<double> dx1;       // delta x (1D array), coordinate unit
        BootesArray<double> dx2;
        BootesArray<double> dx3;
        BootesArray<double> dx1p;      // delta x (1D array), physical unit
        BootesArray<double> dx2p;
        BootesArray<double> dx3p;

        // for spherical polar coordinate, since dxi is not sufficient.
        BootesArray<double> vol;          // volume of each cell
        BootesArray<double> f1a;          // face size in x1 (3D array), size (Nx3, Nx2, Nx1 + 1), Nxi = nxi + 2 * ngi
        BootesArray<double> f2a;          //
        BootesArray<double> f3a;          //
        BootesArray<double> one_orgeo;    // = 0.5 * (rp^2 - rm^2) / (1/3 * rp^3 - rm^3) used in geometry terms
        BootesArray<double> rV;           // = (rm + rp) * 1/3 * (rp^3 - rm^3)
        BootesArray<double> geo_cot;         // = (sin(tp) - sin(tm)) / abs(cos(tp) - cos(tm))
        BootesArray<double> geo_sm;         // sin(tm)
        BootesArray<double> geo_sp;         // sin(tp)
        BootesArray<double> rsq;            // r^2

        double minx1, maxx1, minx2, maxx2, minx3, maxx3, ratio_dim1;
        int x1s, x2s, x3s;                     // start index of active domain
        int x1l, x2l, x3l;                     // end index of active domain
        int nx1, nx2, nx3;                     // number of active zones in each direction
        int ng1, ng2, ng3;                     // number of ghost zones in each direction, implement for 2D and 1D simulation

        double hydro_gamma;
        double vth_coeff;
        #ifdef ENABLE_TEMPERATURE_PROTECTION
        double minTemp;
        #endif // ENABLE_TEMPERATURE_PROTECTION
        #ifdef DENSITY_PROTECTION
        double minDensity;
        #endif

        /** cons **/
        BootesArray<double> cons;           // 4D (5, z, y, x)

        /** prim **/
        BootesArray<double> prim;            // 4D (4, z, y, x)

        /** multi-fluid for dust **/
        int NUMSPECIES;
        double rhodm;           // material density of dust grain. (1D array)
        BootesArray<double> GrainEdgeList;              // edge of dust grains. (1D array, size NUMSPECIES + 1)
        BootesArray<double> GrainSizeList;              // size of dust grains. (1D array)
        BootesArray<double> GrainMassList;              // mass of dust grains. (1D array)
        BootesArray<double> GrainSizeTimesGrainDensity; // rhodm * s
        BootesArray<double> dcons;
        BootesArray<double> dprim;

        /** grav **/
        #if defined (ENABLE_GRAVITY)
            gravity *grav = new gravity;
        #endif
        /** viscosity **/
        #ifdef ENABLE_VISCOSITY
            BootesArray<double> nu_vis;
        #endif // ENABLE_VISCOSITY

        /** setup grid functions **/
        void SetupCartesian(int dimension,
                            double x1min, double x1max, int numx1, int ngh1,
                            double x2min, double x2max, int numx2, int ngh2,
                            double x3min, double x3max, int numx3, int ngh3);
        void SetupSphericalPolar(int dimension,
                                 double x1min, double x1max, int numx1, double ratio1, int ngh1,
                                 double x2min, double x2max, int numx2,                int ngh2,
                                 double x3min, double x3max, int numx3,                int ngh3);
        #if defined (ENABLE_DUSTFLUID)
            void setupDustFluidMesh(int NS);
        #endif

        /** user-defined miscellous quantities **/
        BootesArray<double> UserScalers;
};


#endif
