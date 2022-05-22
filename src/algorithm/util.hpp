#ifndef UTIL_HPP_
#define UTIL_HPP_
#include <cmath>
#include "BootesArray.hpp"

using namespace std;

BootesArray<double> linspace(double vmin, double vmax, int n, bool endpoint){
    BootesArray<double> out;
    out.NewBootesArray(n);
    if (endpoint){
        double dx = (vmax - vmin) / ((double) n - 1.);
        for (int ii = 0; ii < n; ii ++){
            out(ii) = vmin + dx * ii;
        }
    }
    else {
        double dx = (vmax - vmin) / ((double) n);
        for (int ii = 0; ii < n; ii ++){
            out(ii) = vmin + dx * ii;
        }
    }
    return out;
}


BootesArray<double> logspace(double vmin, double vmax, int n, bool endpoint){
    BootesArray<double> out;
    out.NewBootesArray(n);
    if (endpoint){
        double dx = (vmax - vmin) / ((double) n - 1.);
        for (int ii = 0; ii < n; ii ++){
            out(ii) = pow(10, vmin + dx * ii);
        }
    }
    else {
        double dx = (vmax - vmin) / ((double) n);
        for (int ii = 0; ii < n; ii ++){
            out(ii) = pow(10, vmin + dx * ii);
        }
    }
    return out;
}


void cross_prod(double a[3], double b[3], double a_cross_b[3]){
    // Note: cross-product is only properly defined for vectors in R^3 and R^7
    // (see https://www.jstor.org/stable/2323537?seq=2#metadata_info_tab_contents for proof)
    a_cross_b[0] = a[1] * b[2] - a[2] * b[1];
    a_cross_b[1] = a[2] * b[0] - a[0] * b[2];
    a_cross_b[2] = a[0] * b[1] - a[1] * b[0];
}


int searchsorted(double x, BootesArray<double> &xlist){
    /*for (int ii = 0; ii < xlist.shape()[0]; ii++){
        cout << xlist(ii) << '\t';
    }
    cout << endl;*/
    int min_ind = 0;
    int max_ind = xlist.shape()[0] - 1;
    int mid_ind = (min_ind + max_ind) / 2;
    if (xlist.shape()[0] > 2) {
        if (x < xlist(0)){ return 0; }
        if (xlist(0) == x){ return 0; }
        if ((xlist(0) < x) && (x <= xlist(1))){ return 1; }
        if (x >= xlist(max_ind)){return max_ind;}
    }
    while (mid_ind != xlist.shape()[0] - 1 && mid_ind != 0) {
        // cout << xlist(min_ind) << '\t' << xlist(mid_ind) << '\t' << x << '\t' << xlist(mid_ind + 1) << '\t' << xlist(max_ind) << '\t' << min_ind << '\t' << mid_ind << '\t' << max_ind << '\n';
        if (xlist(mid_ind) == x){
            return mid_ind;
        }
        else if ((xlist(mid_ind) < x) && (x <= xlist(mid_ind + 1))){
            return mid_ind + 1;
        }
        else if((xlist(mid_ind) <= x) && (xlist(mid_ind + 1) < x)){
            min_ind = mid_ind;
            mid_ind = (min_ind + max_ind) / 2;
        }
        else if((xlist(mid_ind) > x) && (xlist(mid_ind + 1) >= x)){
            max_ind = mid_ind + 1;
            mid_ind = (min_ind + max_ind) / 2;
        }
        /*else if(xlist(mid_ind + 1) == x){
            return mid_ind + 1;
        }
        else if(xlist(mid_ind) == x){
            return mid_ind;
        }*/
        if (min_ind + 1 >= max_ind){
            cout << ("searchsorted failed") << endl;
            throw ;
        }
    }
    return -1;
}


string choosenumber(int num)
{
    if (num < 10) {
        return "0000" + to_string(num);
    }
    else if (num < 100) {
        return "000" + to_string(num);
    }
    else if (num < 1000) {
        return "00" + to_string(num);
    }
    else if (num < 10000) {
        return "0" + to_string(num);
    }
    else if (num < 100000) {
        return to_string(num);
    }
    else {
        throw "1";
    }
}


void generate_uniform_grav(double unif_grav_acc, unsigned int GridSize[3], BootesArray<double> &Phi,
                           BootesArray<double> &xcoord, BootesArray<double> &ycoord, BootesArray<double> &zcoord){
    double Phi_calc = 0;
    Phi.NewBootesArray(GridSize[2], GridSize[1], GridSize[0]);
    for (int kk = 0; kk < GridSize[2]; kk ++){
        for (int jj = 0; jj < GridSize[1]; jj ++){
            for (int ii = 0; ii < GridSize[0]; ii ++){
                Phi(kk, jj, ii) = Phi_calc;
            }
        }
        if (kk != GridSize[2] - 1){
            Phi_calc += unif_grav_acc * (zcoord(kk + 1) - zcoord(kk));
        }
        else{
            Phi_calc += unif_grav_acc * (zcoord(kk) - zcoord(kk - 1));
        }
    }
}



#endif
