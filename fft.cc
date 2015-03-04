#include <cmath>
#include <algorithm>
#include <functional>
#include "fft.h"
#include <iostream>


const double PI = acos(-1.0);

// performs forward FFT on an input array (one-dimensional) of M
// complex samples.  It is assumed that the input array has size M
// being a power of 2.  The transform is done in-place.
//
// buffer is a temporary array of size M
// start is the start of the array portion to process
// stride is the distance between successive elements in the subarray
//   to process.
void FFT(Array f, int M, Array buffer, int start, int stride)
{
  if (M == 1) {
    // base case, nothing to do
    return;
  }

  // recursively compute the odd and even part
  int K = M/2;

  FFT(f, K, buffer, start, 2*stride);                 // even part
  FFT(f, K, buffer, start + stride, 2*stride);        // odd part

  Coeff omega = exp(Coeff(0, -2*PI/M));
  int ei = start, oi = start + stride;
  Coeff power = 1;
  for (int u = 0; u < K; u++) {
    buffer[u] = f[ei] + f[oi] * power;
    buffer[u+K] = f[ei] - f[oi] * power;
    ei += 2*stride;
    oi += 2*stride;
    power *= omega;
  }

  int i = start;
  for (int u = 0; u < M; u++) {
    f[i] = buffer[u];
    i += stride;
  }

}

// performs inverse FFT on an input array (one-dimensional) of M
// complex samples.  It is assumed that the input array has size M
// being a power of 2.  The transform is done in-place.
//
// buffer is a temporary array of size M
// start is the start of the array portion to process
// stride is the distance between successive elements in the subarray
//   to process.
void inverseFFT(Array f, int M, Array buffer, int start, int stride)
{
  int index = start;
  for (int i = 0; i < M; i++) {
    f[index] = conj(f[index]);
    index += stride;
  }
  FFT(f, M, buffer, start, stride);

  index = start;
  for (int i = 0; i < M; i++) {
    f[index] = conj(f[index])/Coeff(M);
    index += stride;
  }
}

