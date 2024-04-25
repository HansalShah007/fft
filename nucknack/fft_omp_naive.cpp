#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <chrono>
#include <thread>
#include <vector>
#include <complex>

using namespace std;
using cd = complex<double>;

//  g++ -o omp_fft fft_omp.cpp -fopenmp
/*
export OMP_NUM_THREADS=32
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
*/

const double PI = acos(-1);


// Function to find the nearest power of 2
int nearest_pow_2(int n) {
    int power = 1;
    while (power < n) power <<= 1;
    return power;
}

// Function for reversing the bits of an integer
int reverse(int num, int lg_n) {
	int res = 0;
	// TODO why go all the way? halfway seems enough
	for (int i = 0; i < lg_n; i++) {
			if (num & (1 << i))
					res |= 1 << (lg_n - 1 - i);
	}
	return res;
}

vector<cd>& fft(vector<cd> & a, bool invert) {
    /*#pragma omp parallel
    {
		int id = omp_get_thread_num();
		#pragma omp critical
        cout << id << endl;
    }*/

	// Finding the number of steps of butterfly FFTs required
	int n = a.size();
	int lg_n = 0;
	while ((1 << lg_n) < n)
			lg_n++;
	
	// Padding the array with 0s to make it the size of power of 2
	if(n!=(1<<lg_n)){ 
			a.resize(1<<lg_n, 0);
			n = 1 << lg_n;
	}

	// Applying bit reversal on the original array (inplace)
	for (int i = 0; i < n; i++) {
			int reversed_bits = reverse(i, lg_n);
			if (i < reversed_bits)
					swap(a[i], a[reversed_bits]);
	}

	// Iterating through pairs of number of different lengths, starting from 2
	for (int len = 2; len <= n; len <<= 1) {

			// Configuring the principal root of unity that will be used for computing the FFT values for this particular length of pairs
			double ang = 2 * PI / len * (invert ? -1 : 1);
			cd wlen(cos(ang), sin(ang));

			// Iterating through all the possible, non-overlapping pairs for this length
			for (int i = 0; i < n; i += len) {

					// Computing the FFT values for each pair inplace
					cd w(1);
					for (int j = 0; j < len / 2; j++) {
							cd u = a[i+j];
							cd v = a[i+j+len/2] * w;
							a[i+j] = u + v;
							a[i+j+len/2] = u - v;
							w *= wlen;
					}

			}
	}

    // Applying the normalization factor if the function was configured for converting inverse FFT
    if (invert) {
        for (cd & x : a)
            x /= n;
    }

    return a;
}

int main() {
	// Clocking and sanity checks for thread counts
	auto start = std::chrono::steady_clock::now();
    int shared_variable = 0;
    #pragma omp parallel shared(shared_variable)
    {
        // Increment the shared variable with a lock
        #pragma omp atomic
        shared_variable++;
    }

	#pragma omp barrier

	#pragma omp master
	{
		// Run fft and sanity check
		vector<std::complex<double>> complexVector;

		// Add some complex values to the vector
		int range = 1000;
		for (int i = 0; i < range; i++) {
			complexVector.push_back(complex<double>(i, 0));
		}

		vector<cd> complexVectorCopy = complexVector;
		
		vector<complex<double>> sol = fft(fft(complexVector, 0), 1);
		auto end = std::chrono::steady_clock::now();

		double diff = 0.;
		for (int i = 0; i < range; ++i) {
			diff += pow(1.0*complexVectorCopy[i].real() - sol[i].real(), 2);
		}

		cout << "diff: " << diff << endl;

		// Clocking and sanity check
		// cout << "Threads " << shared_variable << " for " << range << endl;
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
		std::cout << "Elapsed time: " << 1.0*duration/1000000 << " sec" << std::endl;
	}
    return 0;
}
