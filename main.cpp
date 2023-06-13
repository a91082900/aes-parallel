#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <mpi.h>
#include "aes.h"

#define DEBUG

using namespace std;

#ifdef MPI
int rk, world_size;
#endif

int main(int argc, char* argv[]) {
    int size = 200;
    if(argc >= 2) {
        size = atoi(argv[1]);
    }

    #ifdef MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if(provided != MPI_THREAD_MULTIPLE) {
        printf("error\n");
        return 0;
    }
    // Get world size and rk
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    #endif

    srand(time(NULL));

    unsigned char key[16] = {
        0xab, 0x9a, 0x8b, 0x7a, 0x6b, 0x5a, 0x4b, 0x3a,
        0x2b, 0x1a, 0x0b, 0xfa, 0xea, 0xda, 0xca, 0xba
        // 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
        // 'i', 'j', 'k','l', 'm', 'n', 'o', 'p'
    };

    chrono::time_point<chrono::steady_clock> t0 = chrono::steady_clock::now();

    unsigned char *in = new unsigned char[size + 16];
    unsigned char *out = new unsigned char[size + 16];
    unsigned char *in_b = new unsigned char[size + 16];
    unsigned char iv[16];

    #ifdef MPI
    if(rk == 0) {
    #endif
        for(int i = 0; i < size; i++) {
            in[i] = rand() % 256;
            in_b[i] = in[i];
        }
        for(int i = 0; i < 16; i++) {
            iv[i] = rand() % 256;
        }
    #ifdef MPI
    }
    #endif
    chrono::time_point<chrono::steady_clock> t1 = chrono::steady_clock::now();
    
    int enc_size = size;
    int counter_init = 0;
    #ifdef MPI
    int sz_per_rank = (size / 16 + 1) / world_size + 1;
    enc_size = sz_per_rank * 16;
    if(rk == 0) {
        for(int r = 1; r < world_size-1; r++) {
            MPI_Send(in + r * enc_size, enc_size, MPI_UNSIGNED_CHAR, r, 0, MPI_COMM_WORLD);
        }
        MPI_Send(in + (world_size-1) * enc_size, size - (world_size-1) * enc_size, MPI_UNSIGNED_CHAR, world_size-1, 0, MPI_COMM_WORLD);
    }
    else if(rk != world_size - 1) {
        MPI_Status status;
        MPI_Recv(in, enc_size, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }
    else {
        enc_size = size - rk * sz_per_rank * 16;
        MPI_Status status;
        MPI_Recv(in, enc_size, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }
    MPI_Bcast(iv, 16, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    counter_init = sz_per_rank * rk;
    #endif
    // encryptECB(in, out, key, size);
    // decryptECB(out, in, key, size);

    // encryptCBC(in, out, key, iv, size);
    // decryptCBC(out, in, key, iv, size);
    
    encryptCTR(in, out, key, iv, enc_size, counter_init);
    encryptCTR(out, in, key, iv, enc_size, counter_init);

    #ifdef MPI
    if(rk == 0) {
        int completed = 1, source, recv_size;
        MPI_Status status;
        while(completed < world_size) {
            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            source = status.MPI_SOURCE;
            if(source != world_size - 1) {
                recv_size = enc_size;
            }
            else {
                recv_size = size - (world_size-1) * enc_size;
            }
            MPI_Recv(out + source * enc_size, recv_size, MPI_UNSIGNED_CHAR, source, 0, MPI_COMM_WORLD, &status);
            completed++;
        }
    }
    else {
        MPI_Send(out, enc_size, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
        return 0;
    }
    #endif

    chrono::time_point<chrono::steady_clock> t2 = chrono::steady_clock::now();
    
    #ifdef DEBUG
    for(int i = 0; i < size; i++) {
        // cout << hex << (int)in[i] << " " << (int)in_b[i] << endl;
        assert(in[i] == in_b[i]);
    }
    #endif
    cout << "Time to generate random input: " << dec << chrono::duration_cast<chrono::microseconds>(t1 - t0).count() << " us" << endl;
    cout << "Time to encrypt and decrypt: " << chrono::duration_cast<chrono::microseconds>(t2 - t1).count() << " us" << endl;

    delete[] in;
    delete[] out;
    delete[] in_b;
    
    return 0;
}