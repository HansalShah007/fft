#include <iostream>
#include <vector>
#include <complex>
#include <cmath>


#include <stdlib.h>
#include <math.h>

#include <fftw3.h>
#include <complex>
#include <cmath>

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>

#define N 4


float sum_milliseconds = 0;
using Clock = std::chrono::high_resolution_clock;

std::chrono::time_point<Clock> start_time, stop_time;
std::chrono::time_point<Clock> start_time_2, stop_time_2;


int main() {

  fftw_complex in[N], out[N], in2[N]; 
  fftw_plan p, q;
  int i;

  for (i = 0; i < N; i++) {
    in[i][0] = cos(3 * 2*M_PI*i/N);
    in[i][1] = 0;
  }
  

  start_time = Clock::now();
  
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  fftw_execute(p);
  
  fftw_destroy_plan(p);


  stop_time = Clock::now();



  int diff = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time).count();
  
  printf("runtime fftw:    n: %d, fft ran in  %d  ms ",  N, diff);

  

  // q = fftw_plan_dft_1d(N, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
  // fftw_execute(q);


  // for (i = 0; i < N; i++) {
  //   in2[i][0] *= 1./N;
  //   in2[i][1] *= 1./N;
  // }
  // // for (i = 0; i < N; i++)
  // //   printf(" %3d %+9.5f %+9.5f I vs. %+9.5f %+9.5f I\n",
  // //       i, in[i][0], in[i][1], in2[i][0], in2[i][1]);
  // fftw_destroy_plan(q);

  fftw_cleanup();



  fftw_complex in3[N], out3[N], in4[N]; 
  fftw_plan p2, q2;
  int i2;

  for (i2 = 0; i2 < N; i2++) {
    in[i2][0] = cos(3 * 2*M_PI*i/N);
    in[i2][1] = 0;
  }
  

  p2 = fftw_plan_dft_1d(N, in3, out3, FFTW_FORWARD, FFTW_ESTIMATE);
  
  start_time_2 = Clock::now();
  fftw_execute(p2);
  stop_time_2 = Clock::now();

  fftw_destroy_plan(p2);

  int diff_2 = std::chrono::duration_cast<std::chrono::microseconds>(stop_time_2 - start_time_2).count();
  
  printf("\n exe only:  %d  ms ",  diff_2);

  fftw_cleanup();
  return 0;


}
