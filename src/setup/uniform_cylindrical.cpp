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
    double Eint;


    void setup(mesh &m, input_file &finput){
        init_unifdensity = 1;
        init_unifinternal = 1;
        //kT_mu_up =  ;
        //kT_mu_low= 0;
        Eint = 1;
        // If you want to read these from the input file at run time, you can use the following two lines
        // init_unifdensity = finput.getDouble("InitUnifDensity");
        // init_unifinternal = finput.getDouble("InitUnifInternal");

        //double temp_slope = (kT_mu_up - kT_mu_low) / (m.x1l - m.x1s);
        //double temp_const = kT_mu_low - m.x3s * temp_slope;
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
     for (int specIND = 0; specIND < m.NUMSPECIES; specIND++){
        for (int kk = m.x3s; kk < m.x3l; kk++){
            for (int jj = m.x2s; jj < m.x2l; jj++){
                for (int ii = m.x1s; ii < m.x1l; ii++){
                    m.dcons(specIND, IDN, kk, jj, ii) = 1.; // (double) specIND + 1.;
                    m.dcons(specIND, IM1, kk, jj, ii) = 0.0;
                    m.dcons(specIND, IM2, kk, jj, ii) = 0.0;
                    m.dcons(specIND, IM3, kk, jj, ii) = 0.0; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                    m.dprim(specIND, IDN, kk, jj, ii) = 1.; // (double) specIND + 1.;
                    m.dprim(specIND, IM1, kk, jj, ii) = 0.0;
                    m.dprim(specIND, IM2, kk, jj, ii) = 0.0;
                    m.dprim(specIND, IM3, kk, jj, ii) = 0.0; // sqrt(G * 10 / pow(100, 3)) * m.x1v(ii);
                }
            }
        }
    }   
    }

    void setup_dust(mesh &m, input_file &finput){
        double smin    = finput.getDouble("smin");
        double smax    = finput.getDouble("smax");
        double rhodm   = finput.getDouble("srho");
        int ns         = finput.getInt("num_species");
        m.NUMSPECIES = ns;
        m.rhodm = rhodm;

        m.GrainEdgeList = logspace(log10(smin), log10(smax), ns + 1, true);
        m.GrainSizeList.NewBootesArray(ns);
        for (int specIND = 0; specIND < m.NUMSPECIES; specIND ++){
            double s2 = m.GrainEdgeList(specIND + 1);
            double s1 = m.GrainEdgeList(specIND);
            m.GrainSizeList(specIND) = pow((pow(s2, 4) - pow(s1, 4)) / (4 * (s2 - s1)), 1./3.);
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
