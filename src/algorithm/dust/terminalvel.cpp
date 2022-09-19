#include "../../defs.hpp"
#include "../mesh/mesh.hpp"
#include "terminalvel.hpp"


#ifdef CARTESIAN_COORD
void dust_terminalvelocityapprixmation_xyz(double &vg1, double &vg2, double &vg3,
                                           double &g1,  double &g2,  double &g3,
                                           double &rhod, double &ts,
                                           double &pd1, double &pd2, double &pd3){
    pd1 = rhod * vg1 + g1 * ts;
    pd2 = rhod * vg2 + g2 * ts;
    pd3 = rhod * vg3 + g3 * ts;
    }
#endif // COORDINATE

#ifdef SPHERICAL_POLAR_COORD
void dust_terminalvelocityapprixmation_rtp(double &vg1, double &vg2, double &vg3,
                                           double &g1,  double &g2,  double &g3,
                                           double &rhod, double &ts,  double &r, double &cottheta,
                                           double &pd1, double &pd2, double &pd3){
    pd1 = rhod * vg1 + (g1 + rhod * (vg2 * vg2 + vg3 * vg3) / r) * ts;
    pd2 = rhod * vg2 + (g2 - rhod * (vg1 * vg2 - vg3 * vg3 * cottheta) / r) * ts;
    pd3 = rhod * vg3 + (g3 - rhod * (vg1 * vg3 + vg2 * vg3 * cottheta) / r) * ts;
    }
#endif // COORDINATE

