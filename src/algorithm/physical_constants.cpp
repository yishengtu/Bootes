#include "physical_constants.hpp"
#include "inoutput/input.hpp"
#include <cmath>


void PhysicalConst::setup_physical_constants(input_file &finput){
    length_scale = finput.getDouble("length_scale");
    time_scale   = finput.getDouble("time_scale");
    mass_scale   = finput.getDouble("mass_scale");

    kb    = cgs_kb / (length_scale * length_scale * mass_scale / (time_scale * time_scale));
    mH    = cgs_mH / mass_scale;
    h     = cgs_h / (length_scale * length_scale * mass_scale / time_scale);
    sigma = cgs_sigma / (mass_scale * pow(length_scale, 4) / time_scale);
    c     = cgs_c / (length_scale / time_scale);
    G     = cgs_G / (pow(length_scale, 3) / (mass_scale * pow(time_scale, 2)));
    year  = cgs_year / time_scale;
    AU    = cgs_AU / length_scale;
    m_sun = cgs_m_sun / mass_scale;
    l_sun = cgs_l_sun / (mass_scale * pow(length_scale, 2) / pow(time_scale, 3));  // 1 solar luminosity
    a     = 4. * sigma / c;
}


