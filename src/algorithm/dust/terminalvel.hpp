#ifndef TERMINALVEL_HPP_
#define TERMINALVEL_HPP_


#ifdef CARTESIAN_COORD
void dust_terminalvelocityapprixmation_xyz(double &vg1, double &vg2, double &vg3,
                                           double &g1,  double &g2,  double &g3,
                                           double &rhod, double &ts,
                                           double &pd1, double &pd2, double &pd3);
#endif // COORDINATE

#ifdef SPHERICAL_POLAR_COORD
void dust_terminalvelocityapprixmation_rtp(double &vg1, double &vg2, double &vg3,
                                           double &g1,  double &g2,  double &g3,
                                           double &rhod, double &ts,  double &r, double &cottheta,
                                           double &pd1, double &pd2, double &pd3);
#endif // COORDINATE


#endif // TERMINALVEL_HPP_



