#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>


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


int main(int argc, char** argv) {
    upcxx::init();
    BUtil::print("HELLO\n");
    if (upcxx::rank_n() > 1) {
	// TODO
    }


    vector<cd> a{ 1, 2, 3, 4 };
    vector<cd> A(4);
    fft(a, 0);
    fft(a, 1);
    for (int i = 0; i < 4; ++i)
        cout << a[i] << "\n";

    /*
    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    std::string kmer_fname = std::string(argv[1]);
    std::string run_type = "";

    if (argc >= 3) {
        run_type = std::string(argv[2]);
    }

    std::string test_prefix = "test";
    if (run_type == "test" && argc >= 4) {
        test_prefix = std::string(argv[3]);
    }

    int ks = kmer_size(kmer_fname);

    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                 "-mers, while this binary is compiled for " +
                                 std::to_string(KMER_LEN) +
                                 "-mers.  Modify packing.hpp and recompile.");
    }

    size_t n_kmers = line_count(kmer_fname);

    // Load factor of 0.5
    size_t hash_table_size = n_kmers * (1.0 / 0.5);
    HashMap hashmap(hash_table_size);

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %d for %d kmers.\n", hash_table_size,
                     n_kmers);
    }

    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());

    if (run_type == "verbose") {
        BUtil::print("Finished reading kmers.\n");
    }

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;

    for (auto& kmer : kmers) {
        bool success = hashmap.insert(kmer);
        if (!success) {
            throw std::runtime_error("Error: HashMap is full!");
        }

        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }
    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();

    double insert_time = std::chrono::duration<double>(end_insert - start).count();
    if (run_type != "test") {
        BUtil::print("Finished inserting in %lf\n", insert_time);
    }
    upcxx::barrier();

    auto start_read = std::chrono::high_resolution_clock::now();

    std::list<std::list<kmer_pair>> contigs;
    for (const auto& start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair kmer;
            bool success = hashmap.find(contig.back().next_kmer(), kmer);
            if (!success) {
                throw std::runtime_error("Error: k-mer not found in hashmap.");
            }
            contig.push_back(kmer);
        }
        contigs.push_back(contig);
    }

    auto end_read = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> read = end_read - start_read;
    std::chrono::duration<double> insert = end_insert - start;
    std::chrono::duration<double> total = end - start;

    int numKmers = std::accumulate(
        contigs.begin(), contigs.end(), 0,
        [](int sum, const std::list<kmer_pair>& contig) { return sum + contig.size(); });

    if (run_type != "test") {
        BUtil::print("Assembled in %lf total\n", total.count());
    }

    if (run_type == "verbose") {
        printf("Rank %d reconstructed %d contigs with %d nodes from %d start nodes."
               " (%lf read, %lf insert, %lf total)\n",
               upcxx::rank_me(), contigs.size(), numKmers, start_nodes.size(), read.count(),
               insert.count(), total.count());
    }

    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        for (const auto& contig : contigs) {
            fout << extract_contig(contig) << std::endl;
        }
        fout.close();
    }
    */
    upcxx::finalize();
    return 0;
}