#include <iostream>
#include <string>
#include "aes.h"

using namespace std;

int main() {
    // string plain = "Hello World!";
    unsigned char key[16] = {
        0xab, 0x9a, 0x8b, 0x7a, 0x6b, 0x5a, 0x4b, 0x3a,
        0x2b, 0x1a, 0x0b, 0xfa, 0xea, 0xda, 0xca, 0xba
        // 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
        // 'i', 'j', 'k','l', 'm', 'n', 'o', 'p'
    };
    // char pl[] = "Hello World!8787"; // 128-bit plain text for single block encryption
    
    // unsigned char state[16];
    // for(int i = 0; i < 16; i++) {
    //     state[i] = pl[i];
    // }
    // encryptBlock(state, key);

    // decryptBlock(state, key);

    string pl = "How are you? I am fine. Thank you!";
    int size = pl.size();
    unsigned char in[size];
    unsigned char out[size + 16];
    for(int i = 0; i < size; i++) {
        in[i] = pl[i];
    }
    encryptECB(in, out, key, size);
    decryptECB(out, in, key, size);
    
    return 0;
}