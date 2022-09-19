#ifndef HLLE_HPP_
#define HLLE_HPP_

#pragma acc routine worker
void hlle( double *valsL,
          double *valsR,
          double *fluxs,
          int IMP,
          double &gamma);

#endif // HLLE_HPP_
