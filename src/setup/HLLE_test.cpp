#include <cstdlib>
#include <cmath>

#include "../algorithm/mesh/mesh.hpp"
#include "../algorithm/BootesArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include "../algorithm/inoutput/input.hpp"
#include "../algorithm/gravity/gravity.hpp"
#include "../algorithm/index_def.hpp"

void setup(mesh &m, input_file &finput){
    double x0 = 0.0;
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                if (m.x1v(ii) < x0){
                    m.cons(IDN, kk, jj, ii) = 1.0;
                    m.cons(IM1, kk, jj, ii) = 0;
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    m.cons(IEN, kk, jj, ii) = 2.;
                }
                else{
                    m.cons(IDN, kk, jj, ii) = 1.0;
                    m.cons(IM1, kk, jj, ii) = 0.;
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    m.cons(IEN, kk, jj, ii) = 1.;
                }
            }
        }
    }
}


void work_after_loop(mesh &m, double &dt){
    ;
}

void setup_dust(mesh &m, input_file &finput){
    ;
}
