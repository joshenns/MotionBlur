#ifndef FFT_H
#define FFT_H

#include <complex>

using std::complex;

// Coefficients are complex numbers
typedef complex<double> Coeff;
typedef Coeff *Array;

// performs forward FFT on an input array (one-dimensional) of M
// complex samples.  It is assumed that the input array has size M
// being a power of 2.  The transform is done in-place.
//	
// buffer is a temporary array of size M
// start is the start of the array portion to process
// stride is the distance between successive elements in the subarray
//   to process.
void FFT(Array f, int M, Array buffer, int start = 0, int stride = 1);

// performs inverse FFT on an input array (one-dimensional) of M
// complex samples.  It is assumed that the input array has size M
// being a power of 2.  The transform is done in-place.
//
// buffer is a temporary array of size M
// start is the start of the array portion to process
// stride is the distance between successive elements in the subarray
//   to process.
void inverseFFT(Array f, int M, Array buffer, int start = 0, int stride = 1);

#endif
