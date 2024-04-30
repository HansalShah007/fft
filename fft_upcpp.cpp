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

#include <chrono>

#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <vector>
#include <complex>
#include <random>

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



float sum_milliseconds = 0;
using Clock = std::chrono::high_resolution_clock;

std::chrono::time_point<Clock> start_time, stop_time;

// Define time point variables for start and stop times
//std::chrono::time_point<Clock> start_time, stop_time;


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



vector<cd>& iterative_fft(vector<cd> & a, bool invert) {

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

    // printf("\n\n\n\n\n\n \n\n");
    // if(upcxx::rank_me() == 0){
        
    //     for (int i = 0; i < a.size(); i++){
    //         printf(" %f \n",  a[i] );
    //     }
    // }

    return a;
}
/*
vector<cd>& fft_old(vector<cd> & a, bool invert) {
    
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
            exit(11);
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
            printf("%d %d  %f\n", upcxx::rank_me(), upcxx::rank_me()*nr_of_elements_per_thread + i, a[upcxx::rank_me()*nr_of_elements_per_thread + i]);
            
            local_ptr[i] = a[upcxx::rank_me()*nr_of_elements_per_thread + i];

            //TODO? use rput? if so: fix wait so it only waits after all rputs were started
            //upcxx::rput(a[0], local_a_ptr[0]).wait(); 
        }

        printf("%f %f \n", local_ptr[0], local_ptr[1]);


        //int nr_of_steps_total = lg_n;

        //iterate through first half of steps 
        for (int s = 1; s <= lg_n/2; s++){
            BUtil::print("rank me %d   step %d \n", upcxx::rank_me(), s);
        }

    return a;
}*/

void sub_fft(vector<cd> & a){
    int n = a.size();
    for (int len = 2; len <= n; len <<= 1) {

        // Configuring the principal root of unity that will be used for computing the FFT values for this particular length of pairs
        double ang = 2 * PI / len * (1);
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
}

void sub_fft_2(vector<cd> & a, int rank, int global_n){
    int n = a.size();

    int step_0 = log2(global_n)/2 ;


    for (int len = 2; len <= n; len <<= 1) {
        int subscript = pow(2, (step_0 + log2(len)));

        // Configuring the principal root of unity that will be used for computing the FFT values for this particular length of pairs
        double ang = 2 * PI / subscript * (1);
        cd wlen_1(cos(ang), sin(ang));
        cd wlen = pow(wlen_1, n);

        //

        // Iterating through all the possible, non-overlapping pairs for this length
        for (int i = 0; i < n; i += len) {

            // Computingh the FFT values for each pair inplace
            //cd w(1);
            cd w = pow(wlen_1, rank );
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }

        }
    }
}

void fft(vector<cd> & a ) {


    int n = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < n)
        lg_n++;
    

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
        exit(11);
    }
    
    //printf("%d ranks %d current \n", upcxx::rank_n(), upcxx::rank_me());
    


    int start = upcxx::rank_me()*  upcxx::rank_n();
    int end = (upcxx::rank_me()+1 )* upcxx::rank_n();

    
    std::vector<complex<double>> local_a_vector(a.begin() + start , a.begin() +  end);
    upcxx::dist_object< std::vector<complex<double>>> local_a = local_a_vector;

    std::vector<cd> local_ref_a = *local_a; //.local();


    // printf(" \n\n ref a before fft  \n\n");
    // for (int i = 0; i < local_a->size(); i++){
    //     printf("%f \n",  local_ref_a[i] );
    // }

    //int nr_of_steps_total = lg_n;


    sub_fft(local_ref_a);



    local_a = local_ref_a;



    // printf(" \n\n ref a after fft  \n\n");
    // for (int i = 0; i < local_a->size(); i++){
    //     printf("%f \n",  local_ref_a[i] );
    // }



    
    upcxx::barrier();

    std::vector<complex<double>>  local_a_copy = local_ref_a;
    //local_a = local_ref_a; 


    

    for (int i = 0; i < local_a->size(); i++){

        if(i != upcxx::rank_me()){
            int target_rank = i;
            cd element = local_a_copy[i];

            upcxx::rpc(target_rank, [] (cd element, int parent_rank,  upcxx::dist_object< std::vector<complex<double>>> &local_a ) {
                
                (*local_a)[parent_rank] = element;
                
            }, element, upcxx::rank_me(), local_a).wait();

        }
    }


    upcxx::barrier();
    local_ref_a = *local_a;


    // printf(" \n\n ref a after communication  \n\n");
    // for (int i = 0; i < local_a->size(); i++){
    //     printf("%f \n",  local_ref_a[i] );
    // }

    




    sub_fft_2(local_ref_a, upcxx::rank_me(), n);



    // local_a = local_ref_a;


    // printf(" \n\n ref a after fft block 2  \n\n");
    // for (int i = 0; i < local_a->size(); i++){
    //     printf("%f \n",  local_ref_a[i] );
    // }





    
    upcxx::barrier();

    local_a_copy = local_ref_a;
    local_a = local_ref_a; 

    for (int i = 0; i < local_a->size(); i++){

        if(i != upcxx::rank_me()){
            int target_rank = i;
            cd element = local_a_copy[i];

            upcxx::rpc(target_rank, [] (cd element, int parent_rank,  upcxx::dist_object< std::vector<complex<double>>> &local_a ) {
                
                (*local_a)[parent_rank] = element;

            }, element, upcxx::rank_me(), local_a).wait();

        }
    }

    upcxx::barrier();
    local_ref_a = *local_a;


    // printf(" \n\n %d ref a after 2nd communication  \n\n", upcxx::rank_me());
    // for (int i = 0; i < local_a->size(); i++){
    //     printf("%f \n",  local_ref_a[i] );
    // }


    upcxx::dist_object< std::vector<complex<double>>> fft_result_all = a;
    if(0 == upcxx::rank_me()){

        std::vector<cd> fft_result = *fft_result_all; 

        for (int i = 0; i < local_a->size(); i++){
            fft_result[i] = local_ref_a[i];
        }

        fft_result_all = fft_result;
    }

    
    upcxx::barrier();

    if(0 != upcxx::rank_me()){
        for (int i = 0; i < local_a->size(); i++){
            cd element = local_ref_a[i];
            upcxx::rpc(0, [] (cd element, int position,  upcxx::dist_object< std::vector<complex<double>>> &fft_result_all ) {
                
                (*fft_result_all)[position] = element;

            }, element, upcxx::rank_me()*nr_of_elements_per_thread + i, fft_result_all).wait();


        }

    }

    upcxx::barrier();

    


    // if(0 == upcxx::rank_me()){
    //     std::vector<cd> fft_result = *fft_result_all; 
    //     printf(" \n\n %d collect on 0  \n\n", upcxx::rank_me());

    //     for (int i = 0; i < n; i++){
    //         printf("%f \n",  fft_result[i] );
    //     }
    // }

}

double random_double() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1000.0);
    return dis(gen);
}

int main(int argc, char** argv) {
    upcxx::init();

    int vector_size = 16;
    std::vector<cd> a(vector_size);
    upcxx::dist_object<std::vector<std::complex<double>>> random_complex_vector_shared;

    
    upcxx::global_ptr<std::complex<double>> b = upcxx::new_array<std::complex<double>>(vector_size);


    if (upcxx::rank_me() == 0) {
        std::vector<std::complex<double>> temp_vector(vector_size);
        for (int i = 0; i < vector_size; ++i) {
            double real_part = static_cast<double>(std::rand()) / RAND_MAX * 1000.0;
            double imag_part = static_cast<double>(std::rand()) / RAND_MAX * 1000.0;
            temp_vector[i] = std::complex<double>(real_part, imag_part);
        }
        std::copy(temp_vector.begin(), temp_vector.end(), b.local());
    }

    upcxx::barrier();

    upcxx::global_ptr<complex<double>> shared_a_ptr = upcxx::broadcast(b, 0).wait();

    std::vector<std::complex<double>> complex_vector(vector_size);
    for (int i = 0; i < vector_size; ++i) {
        complex_vector[i] = upcxx::rget(shared_a_ptr + i).wait();
    }

    a = complex_vector;




   
    upcxx::barrier();


    

    // printf(" \n\n  a  %d  \n\n", upcxx::rank_me());
    // for (int i = 0; i < a.size(); i++){
    //     printf("%f \n",  a[i] );
    // }

    




    start_time = Clock::now();

    int n = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < n)
        lg_n++;
    
    // Padding the array with 0s to make it the size of power of 2
    if(n!=(1<<lg_n)){ 
        a.resize(1<<lg_n, 0);
        n = 1 << lg_n;
    }

    fft(a);







    stop_time = Clock::now();
    
    int diff = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time).count();

    int sum_to_P0 = upcxx::reduce_one(diff, upcxx::op_fast_add, 0).wait();
    //upcxx::global_ptr<complex<double>> shared_a_ptr = upcxx::broadcast(upcxx::new_array<complex<double>>(n), 0).wait();

    //printf("rank  %d  diff   %d", upcxx::rank_me(), diff);

    if (upcxx::rank_me() == 0){

        int avg_recv = 0;
        int rank_n = upcxx::rank_n();

        printf("ranks: %d, n: %d, fft ran in  %d  ms ", upcxx::rank_n(), n, sum_to_P0/upcxx::rank_n());

    }



    //vector<cd> b{ 1, 2, 3, 4,6,7,8,9,10, 11,12,13,14,15,16 };
    //iterative_fft(b, 0);


    
    upcxx::finalize();
}