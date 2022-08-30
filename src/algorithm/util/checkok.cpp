#include "checkok.hpp"
#include "../mesh/mesh.hpp"
#include <cmath>


int check_ok(mesh &m){
    int stat = 0;
    for (int varIND = 0; varIND < m.cons.shape()[0]; varIND ++){
        for (int kk = 0; kk < m.cons.shape()[1]; kk ++){
            for (int jj = 0; jj < m.cons.shape()[2]; jj ++){
                for (int ii = 0; ii < m.cons.shape()[3]; ii ++){
                    if (std::isnan(m.cons(varIND, kk, jj, ii))){
                        std::cout << "nan in m.cons. IND: (" << varIND << ", " << kk << ", " << jj << ", " << ii << ")" << std::endl;
                        stat = 1;
                    }
                    if (std::isnan(m.prim(varIND, kk, jj, ii))){
                        std::cout << "nan in m.prim. IND: (" << varIND << ", " << kk << ", " << jj << ", " << ii << ")" << std::endl;
                        stat = 1;
                    }
                }
            }
        }
    }
    #ifdef ENABLE_DUSTFLUID
    for (int specIND = 0; specIND < m.dcons.shape()[0]; specIND ++){
        for (int varIND = 0; varIND < m.dcons.shape()[1]; varIND ++){
            for (int kk = 0; kk < m.dcons.shape()[2]; kk ++){
                for (int jj = 0; jj < m.dcons.shape()[3]; jj ++){
                    for (int ii = 0; ii < m.dcons.shape()[4]; ii ++){
                        bool printvals = false;
                        if (std::isnan(m.dcons(specIND, varIND, kk, jj, ii))){
                            std::cout << "nan in m.dcons. IND: (" << specIND << ", " << varIND << ", " << kk << ", " << jj << ", " << ii << ")" << std::endl;
                            stat = 1;
                            printvals = true;
                        }
                        if (std::isnan(m.dprim(specIND, varIND, kk, jj, ii))){
                            std::cout << "nan in m.dprim. IND: (" << specIND << ", " << varIND << ", " << kk << ", " << jj << ", " << ii << ")" << std::endl;
                            stat = 1;
                            printvals = true;
                        }
                        if (printvals){
                            std::cout << "hydro cons: " << m.cons(0, kk, jj, ii) << ", " << m.cons(1, kk, jj, ii) << ", " << m.cons(2, kk, jj, ii)
                                                        << ", " << m.cons(3, kk, jj, ii) << ", " << m.cons(4, kk, jj, ii) << std::endl;
                            std::cout << "hydro prim: " << m.prim(0, kk, jj, ii) << ", " << m.prim(1, kk, jj, ii) << ", " << m.prim(2, kk, jj, ii)
                                                        << ", " << m.prim(3, kk, jj, ii) << ", " << m.prim(4, kk, jj, ii) << std::endl;
                            std::cout << "dust dcons: " << m.dprim(specIND, 0, kk, jj, ii) << ", " << m.dprim(specIND, 1, kk, jj, ii) << ", " << m.dprim(specIND, 2, kk, jj, ii)
                                                        << ", " << m.dprim(specIND, 3, kk, jj, ii) << ", " << m.dprim(specIND, 4, kk, jj, ii) << std::endl;
                        }
                    }
                }
            }
        }
    }
    #endif // ENABLE_DUSTFLUID
    std::cout << std::flush;
    return stat;
}
