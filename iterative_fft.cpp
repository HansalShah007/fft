#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;
using cd = complex<double>;
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
    for (int i = 0; i < lg_n; i++) {
        if (num & (1 << i))
            res |= 1 << (lg_n - 1 - i);
    }
    return res;
}

vector<cd>& fft(vector<cd> & a, bool invert) {

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

            // Computingh the FFT values for each pair inplace
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i+j], v = a[i+j+len/2] * w;
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

PYBIND11_MODULE(iterative_fft, m) {
    m.doc() = "Iterative FFT module implemented in C++";
    m.def("fft", &fft, "A function which computes the FFT and IFFT of a complex number vector");
    m.def("nearest_pow_2", &nearest_pow_2, "A function that returns the next nearest value that is a power of 2");
}