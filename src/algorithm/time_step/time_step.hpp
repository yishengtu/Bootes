#ifndef TIME_STEP_HPP_
#define TIME_STEP_HPP_


class mesh;


double timestep(mesh &m, double &CFL);

#endif // TIME_STEP_HPP_
