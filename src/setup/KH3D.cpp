#include "../algorithm/mesh.hpp"
#include "../algorithm/YishengArray.hpp"
#include "../algorithm/boundary_condition/standard_bc.hpp"
#include <cstdlib>
#include "../algorithm/inoutput/input.hpp"

void setup(mesh &m, input_file &finput){
    // step 1: read in necessary parameters from input file
    float gamma_hydro = 1.4;

    int dim     = finput.getInt("dimension");
    float x1min = finput.getFloat("x1min");
    float x1max = finput.getFloat("x1max");
    int nx1     = finput.getInt("nx1");
    float x2min = finput.getFloat("x2min");
    float x2max = finput.getFloat("x2max");
    int nx2     = finput.getInt("nx2");
    float x3min = finput.getFloat("x3min");
    float x3max = finput.getFloat("x3max");
    int nx3     = finput.getInt("nx3");
    // determine number of ghost zones
    int ng1, ng2, ng3;
    if (dim == 1){ ng1 = 2; ng2 = 0; ng3 = 0; }
    else if (dim == 2){ ng1 = 2; ng2 = 2; ng3 = 0; }
    else if (dim == 3){ ng1 = 2; ng2 = 2; ng3 = 2; }
    else { cout << "dimension not recognized! " << endl << flush; throw 1; }

    // step 2: setup grid
    m.SetupCartesian(dim,
                     x1min, x1max, nx1, ng1,                       // ax1
                     x2min, x2max, nx2, ng2,                       // ax2
                     x3min, x3max, nx3, ng3                        // ax3
                     );
    m.hydro_gamma = gamma_hydro;

    // step 3: setup initial condition
    srand(time(0));
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                if (abs(m.x2v(jj)) < 0.1 && abs(m.x3v(kk)) < 0.1){
                    m.cons(IDN, kk, jj, ii) = 2.0;
                    m.cons(IM1, kk, jj, ii) = - 0.5 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    if (m.x1s + m.nx1 * 2 / 5 - 1 < ii && ii < m.x1s + m.nx1 * 3 / 5){
                        m.cons(IM2, kk, jj, ii) += 0.1 * ((float) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                        m.cons(IM3, kk, jj, ii) += 0.1 * ((float) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    }
                    m.prim(IPN, kk, jj, ii) = 2.5;
                    float vel1 = m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel2 = m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel3 = m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    m.cons(IEN, kk, jj, ii) = ene(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii),
                                                                  vel1,
                                                                  vel2,
                                                                  vel3, m.hydro_gamma);
                }
                else{
                    m.cons(IDN, kk, jj, ii) = 1.0;
                    m.cons(IM1, kk, jj, ii) = 0.5 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    if (m.x1s + m.nx1 * 3 / 7 - 1 < ii && ii < m.x1s + m.nx1 * 4 / 7){
                        m.cons(IM2, kk, jj, ii) += 0.1 * ((float) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                        m.cons(IM3, kk, jj, ii) += 0.1 * ((float) rand()/RAND_MAX  - 0.5) * m.cons(IDN, kk, jj, ii);
                    }
                    m.prim(IPN, kk, jj, ii) = 2.5;
                    float vel1 = m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel2 = m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel3 = m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    m.cons(IEN, kk, jj, ii) = ene(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii),
                                                                  vel1,
                                                                  vel2,
                                                                  vel3, m.hydro_gamma);
                }
            }
        }
    }

    // step 4: fill in primitive variables
    m.cons_to_prim();
    // step 5: fill in boundary values in ghost zones
    apply_boundary_condition(m);
}
