#include "../algorithm/mesh.hpp"
#include "../algorithm/YishengArray.hpp"

void setup(mesh &m){
    float gamma_hydro = 1.4;
    m.SetupCartesian(1,
                     0., 1., 100., 2,             // ax1
                     0., 1., 1., 0,               // ax2
                     0., 1., 1., 0                // ax3
                     );
    m.hydro_gamma = gamma_hydro;

    float x0 = 0.3;
    for (int kk = m.x3s; kk < m.x3l; kk++){
        for (int jj = m.x2s; jj < m.x2l; jj++){
            for (int ii = m.x1s; ii < m.x1l; ii++){
                // 1
                if (m.x1v(ii) < x0){
                    m.cons(IDN, kk, jj, ii) = 1.0;
                    m.cons(IM1, kk, jj, ii) = 0.75 * m.cons(IDN, kk, jj, ii);
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    m.prim(IPN, kk, jj, ii) = 1.0;
                    float vel1 = m.cons(IM1, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel2 = m.cons(IM2, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    float vel3 = m.cons(IM3, kk, jj, ii) / m.cons(IDN, kk, jj, ii);
                    m.cons(IEN, kk, jj, ii) = ene(m.cons(IDN, kk, jj, ii), m.prim(IPN, kk, jj, ii),
                                                                  vel1,
                                                                  vel2,
                                                                  vel3, m.hydro_gamma);
                }
                else{
                    m.cons(IDN, kk, jj, ii) = 0.125;
                    m.cons(IM1, kk, jj, ii) = 0;
                    m.cons(IM2, kk, jj, ii) = 0;
                    m.cons(IM3, kk, jj, ii) = 0;
                    m.prim(IPN, kk, jj, ii) = 0.1;
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
    //for (int jj = m.x2s; jj < m.x2l; jj ++){
    //    m.mom1(m.nghost, jj, m.nghost) = 3;
    //}
    //m.mom1(1, m.x2s + m.nx2 * 1 / 4 + 2, 101) = 10;

    m.cons_to_prim();
    apply_boundary_condition(m);
}
