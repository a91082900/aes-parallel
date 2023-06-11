#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>
#include "aes.h"

#define DEBUG

using namespace std;

int main(int argc, char* argv[]) {
    int size = 200;
    if(argc >= 2) {
        size = atoi(argv[1]);
    }
    // srand(time(NULL));

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
    for(int i = 0; i < size; i++) {
        in[i] = rand() % 256;
        in_b[i] = in[i];
    }
    unsigned char iv[16];
    for(int i = 0; i < 16; i++) {
        iv[i] = rand() % 256;
    }

    chrono::time_point<chrono::steady_clock> t1 = chrono::steady_clock::now();
    cout << "Time to generate random input: " << chrono::duration_cast<chrono::microseconds>(t1 - t0).count() << " us" << endl;

    // encryptECB(in, out, key, size);
    // decryptECB(out, in, key, size);

    // encryptCBC(in, out, key, iv, size);
    // decryptCBC(out, in, key, iv, size);
    
    encryptCTR(in, out, key, iv, size);
    encryptCTR(out, in, key, iv, size);

    chrono::time_point<chrono::steady_clock> t2 = chrono::steady_clock::now();
    
    #ifdef DEBUG
    for(int i = 0; i < size; i++) {
        // cout << hex << (int)in[i] << " " << (int)in_b[i] << endl;
        assert(in[i] == in_b[i]);
    }
    #endif

    cout << "Time to encrypt and decrypt: " << dec << chrono::duration_cast<chrono::microseconds>(t2 - t1).count() << " us" << endl;

    delete[] in;
    delete[] out;
    delete[] in_b;
    
    return 0;
}