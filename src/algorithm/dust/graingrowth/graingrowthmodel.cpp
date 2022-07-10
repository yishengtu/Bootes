#include "graingrowthmodel.hpp"
#include <cmath>

void grain_growth_model_stick(double &s1, double &s2, double &dv, double res[2]){
    res[0] = dv * M_PI * pow((s1 + s2), 2.0);
    res[1] = 0.0;
}
