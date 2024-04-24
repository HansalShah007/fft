#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <chrono>
#include <thread>
#include <vector>
#include <complex>

using namespace std;

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

// Recursive FFT function
vector<complex<double>> fft(vector<complex<double>>& a) {
    int n = a.size();

    if (n == 1)
        return vector<complex<double>>(1, a[0]);

    // Check if the size is a power of 2, if not pad with zeros
    int n2 = nearest_pow_2(n);
    if (n2 > n) {
        a.resize(n2, 0);
        n = n2;
    }

    vector<complex<double>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    auto y_even = fft(even);
    auto y_odd = fft(odd);

    vector<complex<double>> y(n);
    complex<double> omega_n = exp(complex<double>(0, 2 * PI / n));
    complex<double> omega = 1;

    for (int j = 0; j < n / 2; ++j) {
        y[j] = y_even[j] + omega * y_odd[j];
        y[j + n / 2] = y_even[j] - omega * y_odd[j];
        omega *= omega_n;
    }

    return y;
}

// Recursive IFFT function
vector<complex<double>> ifft(vector<complex<double>>& a, int depth) {
    int n = a.size();

    if (n == 1)
        return vector<complex<double>>(1, a[0]);

    int n2 = nearest_pow_2(n);
    if (n2 > n) {
        a.resize(n2, 0);
        n = n2;
    }

    vector<complex<double>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    auto y_even = ifft(even, depth+1);
    auto y_odd = ifft(odd,depth+1);

    vector<complex<double>> y(n);
    complex<double> omega_n = exp(complex<double>(0, -2 * PI / n));
    complex<double> omega = 1;

    for (int j = 0; j < n / 2; ++j) {
        y[j] = y_even[j] + omega * y_odd[j];
        y[j + n / 2] = y_even[j] - omega * y_odd[j];
        omega *= omega_n;
    }

    if (depth==1){
        for (auto& value : y) value /= y.size();
    }

    return y;
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
		int range = 10000000;
		for (int i = 0; i < range; i++) {
			complexVector.push_back(complex<double>(5*sin(i), 0));
		}
		
		vector<complex<double>> returned = fft(complexVector);
		vector<complex<double>> sol = ifft(returned, 1);

		double diff = 0.;
		for (int i = 0; i < range; ++i) {
			diff += pow(complexVector[i].real() - sol[i].real(), 2);
		}

		cout << "diff: " << diff << endl;

		// Clocking and sanity check
		cout << "Threads " << shared_variable << " for " << range << endl;
		auto end = std::chrono::steady_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
		std::cout << "Elapsed time: " << 1.0*duration/1000000 << " sec" << std::endl;
	}
    return 0;
}
