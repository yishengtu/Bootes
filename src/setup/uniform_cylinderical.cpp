#include <cstdlib>
#include <cmath>

#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include "../algorithm/inoutput/input.hpp"
#include "../algorithm/gravity/gravity.hpp"
#include "../algorithm/index_def.hpp"


/** Setup a uniform initial condition **/
namespace {
    double init_unifdensity;
    double init_unifinternal;


    void setup(mesh &m, input_file &finput){
        init_unifdensity = 1;
        init_unifinternal = 1;
        // If you want to read these from the input file at run time, you can use the following two lines
        // init_unifdensity = finput.getDouble("InitUnifDensity");
        // init_unifinternal = finput.getDouble("InitUnifInternal");

        double temp_slope = (kT_mu_up - kT_mu_low) / (m.x1l - m.x1s);
        double temp_const = kT_mu_low - m.x3s * temp_slope;
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    double rho = init_unifdensity;
                    m.cons(IDN, kk, jj, ii) = rho;
                    m.cons(IM1, kk, jj, ii) = 0.0;
                    m.cons(IM2, kk, jj, ii) = 0.0;
                    m.cons(IM3, kk, jj, ii) = 0.0;
                    m.cons(IEN, kk, jj, ii) = Eint;
                }
            }
        }
    }


    void work_after_loop(mesh &m, double &dt){
        ;
    }

    void apply_user_extra_boundary_condition(mesh &m){
        ;
    }


    void calculate_nu_vis(mesh &m){
        ;
    }
}
