#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <memory>


/*
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>*/
#include <complex>
#include <cmath>


using namespace std;
 
typedef complex<double> cd;

const double PI = acos(-1);

#include "butil.hpp"




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
    
        //TODO make sure this isn't computed on every rank for now


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

        //write  into reversed a at a position even if not on rank  - lots of communication if done on separate ranks
        // Applying bit reversal on the original array (inplace)
        for (int i = 0; i < n; i++) {
            int reversed_bits = reverse(i, lg_n);
            if (i < reversed_bits)
                swap(a[i], a[reversed_bits]);
        }



        int nr_of_elements_per_thread = sqrt(n);
        int nr_of_threads = sqrt(n);


        if (nr_of_threads != upcxx::rank_n()){
            exit(1);
        }


        //TODO make a on rank 0 
        // if (upcxx::rank_me() == 0) {
        //    shared_a_ptr = upcxx::new_array<complex<double>>(n); 
        // }

        //TODO fix this, I'm pretty sure rank_me() is always 0
        upcxx::global_ptr<complex<double>> shared_a_ptr = upcxx::broadcast(upcxx::new_array<complex<double>>(n), 0).wait();
        upcxx::global_ptr<complex<double>> local_a_ptr = shared_a_ptr + upcxx::rank_me()*nr_of_elements_per_thread;

        complex<double> *local_ptr = local_a_ptr.local();

        for (int i = 0; i< nr_of_elements_per_thread; i++){
            //TODO this definitely isn't iterating through correctly
            BUtil::print("%d %d  %f\n", upcxx::rank_me(), upcxx::rank_me()*nr_of_elements_per_thread + i, a[upcxx::rank_me()*nr_of_elements_per_thread + i]);
            
            local_ptr[i] = a[upcxx::rank_me()*nr_of_elements_per_thread + i];

            //TODO? use rput? if so: fix wait so it only waits after all rputs were started
            //upcxx::rput(a[0], local_a_ptr[0]).wait(); 
        }

        BUtil::print("%f %f \n", local_ptr[0], local_ptr[1]);


        //int nr_of_steps_total = lg_n;

        //iterate through first half of steps 
        for (int s = 1; s <= lg_n/2; s++){
            BUtil::print("rank me %d   step %d \n", upcxx::rank_me(), s);
        }

    return a;
}


int main(int argc, char** argv) {
    upcxx::init();
    
    BUtil::print("%d ranks %d current \n", upcxx::rank_n(), upcxx::rank_me());
    
    vector<cd> a{ 1, 2, 3, 4 };

    fft(a, 0);
    fft(a, 1);

    for (int i = 0; i < 4; ++i)
        cout << a[i] << "\n";
    
    return 0;
    
    upcxx::finalize();
}