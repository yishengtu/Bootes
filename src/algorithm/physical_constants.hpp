#ifndef PHYSICAL_CONSTANTS_HPP_
#define PHYSICAL_CONSTANTS_HPP_
#include <cmath>

const double length_scale = 1;
const double time_scale = 1;
const double mass_scale = 1;

const double cgs_kb = 1.380649 * pow(10, -16); // Boltzmann constant
const double cgs_mH = 1.6733 * pow(10, -24);   // mass of neutral hydrogen
const double cgs_h = 6.6261 * pow(10, -27);    // plank's constant
const double cgs_sigma = 5.6704 * pow(10, -5); // Stefan-Boltzmann constant
const double cgs_c = 2.99792458 * pow(10, 10);   // light speed
const double cgs_G = 6.674 * pow(10, -8);      // gravitational constant
const double cgs_year = 31553241;  // 31536000;              // 1 year = 365 days * 24 hrs / day * 3600 sec / hrs
const double cgs_AU = 1.496 * pow(10, 13);     // 1AU
const double cgs_m_sun = 1.989 * pow(10, 33);  // solar mass
const double cgs_l_sun = 3.839 * pow(10, 33);  // 1 solar luminosity
const double cgs_a = 7.5646 * pow(10, -15);    // radiation constant

const float hydro_mu = 2.33;
const float hydro_gamma = 5.0/3.0;

const double kb    = cgs_kb / (length_scale * length_scale * mass_scale / (time_scale * time_scale));
const double mH    = cgs_mH / mass_scale;
const double h     = cgs_h / (length_scale * length_scale * mass_scale / time_scale);
const double sigma = cgs_sigma / (mass_scale * pow(length_scale, 4) / time_scale);
const double c     = cgs_c / (length_scale / time_scale);
const double G     = cgs_G / (pow(length_scale, 3) / (mass_scale * pow(time_scale, 2)));
const double year  = cgs_year / time_scale;
const double AU    = cgs_AU / length_scale;
const double m_sun = cgs_m_sun / mass_scale;
const double l_sun = cgs_l_sun / (mass_scale * pow(length_scale, 2) / pow(time_scale, 3));  // 1 solar luminosity
const double a     = 4. * sigma / c;

#endif
