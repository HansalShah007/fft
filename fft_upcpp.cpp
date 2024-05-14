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


float sum_communication = 0;
using Clock = std::chrono::high_resolution_clock;

std::chrono::time_point<Clock> start_time, stop_time;
std::chrono::time_point<Clock> start_time_2, stop_time_2;
std::chrono::time_point<Clock> start_time_3, stop_time_3;
std::chrono::time_point<Clock> start_time_4, stop_time_4;


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


vector<complex<double>> recursive_fft(vector<complex<double>>& a) {
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

    auto y_even = recursive_fft(even);
    auto y_odd = recursive_fft(odd);

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

    printf("\n\n\n\n\n\n \n\n");
    if(upcxx::rank_me() == 0){
        for (int i = 0; i < a.size(); i++){
            printf(" %f \n",  a[i] );
        }
    }


    return a;
}


void sub_fft(upcxx::global_ptr<std::complex<double>> a, int n){
    
    for (int len = 2; len <= n; len <<= 1) {

        // Configuring the principal root of unity that will be used for computing the FFT values for this particular length of pairs
        double ang = 2 * PI / len * (1);
        cd wlen(cos(ang), sin(ang));

        // Iterating through all the possible, non-overlapping pairs for this length
        for (int i = 0; i < n; i += len) {

            // Computingh the FFT values for each pair inplace
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = upcxx::rget(a+(i+j)).wait();
                cd v = upcxx::rget(a+(i+j+len/2)).wait() * w;
                upcxx::rput(u+v, a+(i+j)).wait();
                upcxx::rput(u-v, a+(i+j+len/2)).wait();
                w *= wlen;
            }

        }
    }
}


void sub_fft_2(upcxx::global_ptr<std::complex<double>> a, int rank, int global_n, int n){
   
    int step_0 = log2(global_n)/2 ;

    for (int len = 2; len <= n; len <<= 1) {
        int subscript = pow(2, (step_0 + log2(len)));

        double ang = 2 * PI / subscript * (1);
        cd wlen_1(cos(ang), sin(ang));
        cd wlen = pow(wlen_1, n);

        
        for (int i = 0; i < n; i += len) {

            cd w = pow(wlen_1, rank );
            for (int j = 0; j < len / 2; j++) {

                cd u = upcxx::rget(a+(i+j)).wait();
                cd v = upcxx::rget(a+(i+j+len/2)).wait() * w;
                upcxx::rput(u+v, a+(i+j)).wait();
                upcxx::rput(u-v, a+(i+j+len/2)).wait();
                w *= wlen;


            }

        }
    }
}




void fft_ideal_ranks(vector<cd> & a) {

    // Calculating the number of steps
    int total_elements = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < total_elements)
        lg_n++;

    // Elements per rank
    int total_ranks = upcxx::rank_n();
    int elements_per_rank = total_elements/total_ranks;

    // All the ranks do the bit reversal
    //write  into reversed a at a position even if not on rank  - lots of communication if done on separate ranks
    // Applying bit reversal on the original array (inplace)
    for (int i = 0; i < total_elements; i++) {
        int reversed_bits = reverse(i, lg_n);
        if (i < reversed_bits)
            swap(a[i], a[reversed_bits]);
    }

    // Temporary check to make sure we are operating on the right number of elements and threads
    int nr_of_elements_per_thread = sqrt(total_elements);
    int nr_of_threads = sqrt(total_elements);

    if (nr_of_threads != total_ranks){
        exit(11);
    }

    // Determining the slice of the array that belongs to the rank
    int start = upcxx::rank_me()*(elements_per_rank);
    int end = (upcxx::rank_me()+1)*(elements_per_rank);

    // Initializing the local data and received data distributed objects
    using dobj_data = upcxx::dist_object<upcxx::global_ptr<complex<double>>>;
    dobj_data local_data = dobj_data(upcxx::new_array<complex<double>>(elements_per_rank));
    dobj_data received_data = dobj_data(upcxx::new_array<complex<double>>(elements_per_rank)); 

    std::vector<complex<double>> local_data_vector(a.begin() + start , a.begin() +  end);

    // Get local access to the global pointer
    upcxx::global_ptr<std::complex<double>> local_ref_data = *local_data;

    // Copy from local vector to global memory
    upcxx::rput(local_data_vector.data(), local_ref_data, end - start).wait();

    // First phase FFT
    sub_fft(local_ref_data, elements_per_rank);
    
    upcxx::barrier();

    // Rearranging the data
    std::vector<complex<double>> local_data_copy(elements_per_rank);
    for (int i = 0; i<elements_per_rank; i++){
        local_data_copy[i] = upcxx::rget(local_ref_data+i).wait();
    }

    start_time_3 = Clock::now();
    for (int i = 0; i < elements_per_rank; i++){

        if(i != upcxx::rank_me()){
            int target_rank = i;
            cd element = local_data_copy[i];

            upcxx::rpc(target_rank, [] (cd element, int parent_rank,  upcxx::dist_object<upcxx::global_ptr<complex<double>>> &local_data) {
                upcxx::rput(element, *local_data + parent_rank).wait();
            }, element, upcxx::rank_me(), local_data).wait();

        }
    }
    stop_time_3 = Clock::now();

    upcxx::barrier();

    
    sub_fft_2(local_ref_data, upcxx::rank_me(), total_elements, elements_per_rank);


    upcxx::barrier();

    // Rearranging for getting the correct order of elements
    for (int i = 0; i<elements_per_rank; i++){
        local_data_copy[i] = upcxx::rget(local_ref_data+i).wait();
    }

    start_time_4 = Clock::now();

    for (int i = 0; i < elements_per_rank; i++){

        if(i != upcxx::rank_me()){
            int target_rank = i;
            cd element = local_data_copy[i];

            upcxx::rpc(target_rank, [] (cd element, int parent_rank,  upcxx::dist_object<upcxx::global_ptr<complex<double>>> &local_data ) {                
                upcxx::rput(element, *local_data + parent_rank).wait();
            }, element, upcxx::rank_me(), local_data).wait();

        }
    }
    stop_time_4 = Clock::now();

    upcxx::barrier();

}





void fft(vector<cd> & a) {

    // Calculating the number of steps
    int total_elements = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < total_elements)
        lg_n++;

    // Elements per rank
    int total_ranks = upcxx::rank_n();
    int elements_per_rank = total_elements/total_ranks;

    // All the ranks do the bit reversal
    //write  into reversed a at a position even if not on rank  - lots of communication if done on separate ranks
    // Applying bit reversal on the original array (inplace)
    for (int i = 0; i < total_elements; i++) {
        int reversed_bits = reverse(i, lg_n);
        if (i < reversed_bits)
            swap(a[i], a[reversed_bits]);
    }


    // Determining the slice of the array that belongs to the rank
    int start = upcxx::rank_me()*(elements_per_rank);
    int end = (upcxx::rank_me()+1)*(elements_per_rank);

    // Initializing the local data and received data distributed objects
    using dobj_data = upcxx::dist_object<upcxx::global_ptr<complex<double>>>;
    dobj_data local_data = dobj_data(upcxx::new_array<complex<double>>(elements_per_rank));
    dobj_data received_data = dobj_data(upcxx::new_array<complex<double>>(elements_per_rank)); 

    std::vector<complex<double>> local_data_vector(a.begin() + start , a.begin() +  end);
    
    // Get local access to the global pointer
    upcxx::global_ptr<std::complex<double>> local_ref_data = *local_data;

    // Copy from local vector to global memory
    upcxx::rput(local_data_vector.data(), local_ref_data, end - start).wait();

    // First phase FFT
    sub_fft(local_ref_data, elements_per_rank);
    
    upcxx::barrier();

    // Calculating the number of steps
    int steps = 0;
    while ((1 << steps) < total_ranks)
        steps++;

    int ranks_per_group = 1;

    for(int i=0; i<steps; i++){

        // Logic for figuring out what rank to send the data to
        ranks_per_group*=2;
        int group_index = upcxx::rank_me()/ranks_per_group;
        int group_local_index = upcxx::rank_me()%ranks_per_group;
        int who_to_send;
        if(group_local_index<(ranks_per_group/2)){
            who_to_send = group_local_index + (ranks_per_group/2) + (group_index*ranks_per_group);
        }
        else {
            who_to_send = group_local_index - (ranks_per_group/2) + (group_index*ranks_per_group);
        }

        // Pack the data into a vector
        std::vector<cd> send_data(elements_per_rank);
        for (int j = 0; j < elements_per_rank; j++) {
            send_data[j] = upcxx::rget(*local_data + j).wait();
        }


        // Send the vector using RPC
        upcxx::rpc(who_to_send, [](std::vector<cd> received_chunk, int sender_rank, upcxx::dist_object<upcxx::global_ptr<cd>>& received_data) {
            upcxx::global_ptr<cd> local_received_data = *received_data;
            for (size_t j = 0; j < received_chunk.size(); j++) {
                upcxx::rput(received_chunk[j], local_received_data + j).wait();
            }
        }, send_data, upcxx::rank_me(), received_data).wait();

       
       // Waiting for all the ranks to finish sending data to each other
        upcxx::barrier();

        // Performing the FFT computation for a single step
        int len = ranks_per_group*elements_per_rank;
        double ang = 2 * PI / len * (1);
        cd wlen(cos(ang), sin(ang));
        cd w(1);

        if(group_local_index<(ranks_per_group/2)){
            w = pow(wlen, elements_per_rank*group_local_index);
        }
        else{
            w = pow(wlen, elements_per_rank*(group_local_index-(ranks_per_group/2)));
        }
        
        for (int j = 0; j < elements_per_rank; j++){
            if(group_local_index<(ranks_per_group/2)){
                cd u = upcxx::rget(*local_data+j).wait();
                cd v = upcxx::rget(*received_data+j).wait() * w;
                upcxx::rput(u+v,*local_data+j).wait();
            }
            else{
                cd u = upcxx::rget(*received_data+j).wait();
                cd v = upcxx::rget(*local_data+j).wait() * w;
                upcxx::rput(u-v,*local_data+j).wait();
            }
            w*=wlen;
        }

        // Waiting for all the ranks to finish their FFT calculations
        upcxx::barrier();

    }

    upcxx::barrier();

    // printf("\n\n%d ref a after 2nd communication\n\n", upcxx::rank_me());
    // for (int i = 0; i < elements_per_rank; i++){
    //     printf("%f\n", upcxx::rget(*local_data+i).wait());
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

    int vector_size = pow(2, 24); 
    std::vector<cd> a(vector_size);
    std::vector<cd> c(vector_size);
    upcxx::dist_object<std::vector<std::complex<double>>> random_complex_vector_shared;

    
    upcxx::global_ptr<std::complex<double>> b = upcxx::new_array<std::complex<double>>(vector_size);


    if (upcxx::rank_me() == 0) {
        std::vector<std::complex<double>> temp_vector(vector_size);
        for (int i = 0; i < vector_size; ++i) {
            double real_part = static_cast<double>(std::rand()) / RAND_MAX * 1000.0;
            double imag_part = static_cast<double>(std::rand()) / RAND_MAX * 1000.0;
            temp_vector[i] = std::complex<double>(real_part, imag_part);
            c[i] = std::complex<double>(real_part, imag_part);
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
    //fft_ideal_ranks(a);
    //iterative_fft(a, 0);
    //recursive_fft(a);

    stop_time = Clock::now();
    
    int diff = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time).count();
    int sum_to_P0 = upcxx::reduce_one(diff, upcxx::op_fast_add, 0).wait();
    
    if (upcxx::rank_me() == 0){

        int avg_recv = 0;
        int rank_n = upcxx::rank_n();

        printf("runtime:   ranks: %d, n: %d, fft ran in  %d  µs ", upcxx::rank_n(), vector_size, sum_to_P0/upcxx::rank_n());

    }

    // int communication_diff = std::chrono::duration_cast<std::chrono::microseconds>(stop_time_3 - start_time_3).count() +  std::chrono::duration_cast<std::chrono::microseconds>(stop_time_4 - start_time_4).count();
    // int sum_to_P0_comm = upcxx::reduce_one(communication_diff, upcxx::op_fast_add, 0).wait();
    
    // if (upcxx::rank_me() == 0){

    //     printf("\n communication time :   %d  µs ",  sum_to_P0_comm/upcxx::rank_n());

    // }
    
    upcxx::finalize();
}