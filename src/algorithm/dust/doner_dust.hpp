#ifndef DONER_DUST_HPP_
#define DONER_DUST_HPP_


// in order to call a function, it needs to be compiled on device
// acc routine: create a device routine
// seq: defines the parallization used within this loop
#pragma acc routine seq
void doner_cell_dust( double *valsL,
          double *valsR,
          double *fluxs,
          int IMP,
          double &gamma);

#endif // DONER_DUST_HPP_
