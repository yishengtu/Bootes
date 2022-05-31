#ifndef UTIL_HPP_
#define UTIL_HPP_
#include "../BootesArray.hpp"

using namespace std;

BootesArray<double> linspace(double vmin, double vmax, int n, bool endpoint);

BootesArray<double> logspace(double vmin, double vmax, int n, bool endpoint);

void cross_prod(double a[3], double b[3], double a_cross_b[3]);

int searchsorted(double x, BootesArray<double> &xlist);

string choosenumber(int num);

#endif
