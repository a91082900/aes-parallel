#include "aes.h"

void addRoundKey(unsigned char* state, unsigned char* key) {
    for(int i = 0; i < 16; i++) {
        state[i] ^= key[i];
    }

    #ifdef DEBUG
    std::cout << endl << "---AddRoundKey Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)state[i] << " ";
    }
    std::cout << endl << "---AddRoundKey End---" << endl;
    #endif
}

void subBytes(unsigned char* state) {
    for(int i = 0; i < 16; i++) {
        state[i] = sbox[state[i]];
    }

    #ifdef DEBUG
    std::cout << endl << "---SubBytes Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)state[i] << " ";
    }
    std::cout << endl << "---SubBytes---" << endl;
    #endif
}

void mixColumns(unsigned char* state) {
    unsigned char tmp[16];
    for(int i = 0; i < 16; i++) {
        tmp[i] = state[i];
    }
    state[0] = GF_2[tmp[0]] ^ GF_3[tmp[1]] ^ tmp[2] ^ tmp[3];
    state[1] = tmp[0] ^ GF_2[tmp[1]] ^ GF_3[tmp[2]] ^ tmp[3];
    state[2] = tmp[0] ^ tmp[1] ^ GF_2[tmp[2]] ^ GF_3[tmp[3]];
    state[3] = GF_3[tmp[0]] ^ tmp[1] ^ tmp[2] ^ GF_2[tmp[3]];
    state[4] = GF_2[tmp[4]] ^ GF_3[tmp[5]] ^ tmp[6] ^ tmp[7];
    state[5] = tmp[4] ^ GF_2[tmp[5]] ^ GF_3[tmp[6]] ^ tmp[7];
    state[6] = tmp[4] ^ tmp[5] ^ GF_2[tmp[6]] ^ GF_3[tmp[7]];
    state[7] = GF_3[tmp[4]] ^ tmp[5] ^ tmp[6] ^ GF_2[tmp[7]];
    state[8] = GF_2[tmp[8]] ^ GF_3[tmp[9]] ^ tmp[10] ^ tmp[11];
    state[9] = tmp[8] ^ GF_2[tmp[9]] ^ GF_3[tmp[10]] ^ tmp[11];
    state[10] = tmp[8] ^ tmp[9] ^ GF_2[tmp[10]] ^ GF_3[tmp[11]];
    state[11] = GF_3[tmp[8]] ^ tmp[9] ^ tmp[10] ^ GF_2[tmp[11]];
    state[12] = GF_2[tmp[12]] ^ GF_3[tmp[13]] ^ tmp[14] ^ tmp[15];
    state[13] = tmp[12] ^ GF_2[tmp[13]] ^ GF_3[tmp[14]] ^ tmp[15];
    state[14] = tmp[12] ^ tmp[13] ^ GF_2[tmp[14]] ^ GF_3[tmp[15]];
    state[15] = GF_3[tmp[12]] ^ tmp[13] ^ tmp[14] ^ GF_2[tmp[15]];

    #ifdef DEBUG
    std::cout << endl << "---MixColumn Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)state[i] << " ";
    }
    std::cout << endl << "---MixColumn End---" << endl;
    #endif
}

void shiftRows(unsigned char *state) {
    unsigned char tmp[16];
    for(int i = 0; i < 16; i++) {
        tmp[i] = state[i];
    }

    state[0] = tmp[0];
    state[1] = tmp[5];
    state[2] = tmp[10];
    state[3] = tmp[15];
    state[4] = tmp[4];
    state[5] = tmp[9];
    state[6] = tmp[14];
    state[7] = tmp[3];
    state[8] = tmp[8];
    state[9] = tmp[13];
    state[10] = tmp[2];
    state[11] = tmp[7];
    state[12] = tmp[12];
    state[13] = tmp[1];
    state[14] = tmp[6];
    state[15] = tmp[11];

    #ifdef DEBUG
    std::cout << endl << "---ShiftRow Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)state[i] << " ";
    }
    std::cout << endl << "---ShiftRow End---" << endl;
    #endif
}

void rotWord(unsigned char *word) {
    unsigned char tmp = word[0];
    word[0] = word[1];
    word[1] = word[2];
    word[2] = word[3];
    word[3] = tmp;
}

void subWord(unsigned char *word) {
    #pragma GCC unroll 4
    for(int i = 0; i < 4; i++) {
        word[i] = sbox[word[i]];
    }
}

void keyExpansion(unsigned char* key, unsigned char* w) {
    unsigned char tmp[4];
    
    for(int i = 0; i < 16; i++) {
        w[i] = key[i];
    }
    for(int i = 16; i < 4 * 4 * 11; i += 4) {
        tmp[0] = w[i - 4];
        tmp[1] = w[i - 3];
        tmp[2] = w[i - 2];
        tmp[3] = w[i - 1];

        if(i % 16 == 0) {
            rotWord(tmp);
            subWord(tmp);
            tmp[0] ^= rc[i / 16]; // rcon[1..3] = 0
        }
        #pragma GCC unroll 4
        for(int j = 0; j < 4; j++) {
            w[i + j] = w[i - 16 + j] ^ tmp[j];
        }
    }

    #ifdef DEBUG
    std::cout << endl << "---KeyExpansion Begin---" << endl;
    for(int i = 0; i < 11; i++) {
        for(int j = 0; j < 16; j++) {
            std::cout << std::hex << (int)w[i * 16 + j] << " ";
        }
        std::cout << endl;
    }
    std::cout << "---KeyExpansion End---" << endl;
    #endif
}

void invSubBytes(unsigned char* state) {
    #pragma GCC unroll 16
    for(int i = 0; i < 16; i++) {
        state[i] = inv_sbox[state[i]];
    }

    #ifdef DEBUG
    std::cout << endl << "---InvSubBytes Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)state[i] << " ";
    }
    std::cout << endl << "---InvSubBytes End---" << endl;
    #endif
}

void invShiftRows(unsigned char* state) {
    unsigned char tmp[16];
    for(int i = 0; i < 16; i++) {
        tmp[i] = state[i];
    }

    state[0] = tmp[0];
    state[1] = tmp[13];
    state[2] = tmp[10];
    state[3] = tmp[7];
    state[4] = tmp[4];
    state[5] = tmp[1];
    state[6] = tmp[14];
    state[7] = tmp[11];
    state[8] = tmp[8];
    state[9] = tmp[5];
    state[10] = tmp[2];
    state[11] = tmp[15];
    state[12] = tmp[12];
    state[13] = tmp[9];
    state[14] = tmp[6];
    state[15] = tmp[3];

    #ifdef DEBUG
    std::cout << endl << "---InvShiftRow Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)state[i] << " ";
    }
    std::cout << endl << "---InvShiftRow End---" << endl;
    #endif
}

void invMixColumns(unsigned char* state) {
    unsigned char tmp[16];
    for(int i = 0; i < 16; i++) {
        tmp[i] = state[i];
    }

    #pragma GCC unroll 4
    for(int i = 0; i < 4; i++) {
        state[4*i + 0] = GF_14[tmp[4*i + 0]] ^ GF_11[tmp[4*i + 1]] ^ GF_13[tmp[4*i + 2]] ^ GF_9[tmp[4*i + 3]];
        state[4*i + 1] = GF_9[tmp[4*i + 0]] ^ GF_14[tmp[4*i + 1]] ^ GF_11[tmp[4*i + 2]] ^ GF_13[tmp[4*i + 3]];
        state[4*i + 2] = GF_13[tmp[4*i + 0]] ^ GF_9[tmp[4*i + 1]] ^ GF_14[tmp[4*i + 2]] ^ GF_11[tmp[4*i + 3]];
        state[4*i + 3] = GF_11[tmp[4*i + 0]] ^ GF_13[tmp[4*i + 1]] ^ GF_9[tmp[4*i + 2]] ^ GF_14[tmp[4*i + 3]];
    }

    #ifdef DEBUG
    std::cout << endl << "---InvMixColumn Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)state[i] << " ";
    }
    std::cout << endl << "---InvMixColumn End---" << endl;
    #endif
}

void encryptBlock(unsigned char* in, unsigned char* out, unsigned char* w) {
    // unsigned char w[4 * 4 * 11];
    // keyExpansion(key, w);

    for(int i = 0; i < 16; i++) {
        out[i] = in[i];
    }

    addRoundKey(out, w);
    for(int i = 1; i < 10; i++) {
        subBytes(out);
        shiftRows(out);
        mixColumns(out);
        addRoundKey(out, w + i * 16);
    }
    subBytes(out);
    shiftRows(out);
    addRoundKey(out, w + 10 * 16);

    #ifdef DEBUG
    std::cout << endl << "---EncryptBlock Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)out[i] << " ";
    }
    std::cout << endl << "---EncryptBlock End---" << endl;
    #endif
}

void decryptBlock(unsigned char* in, unsigned char* out, unsigned char* w) {
    // unsigned char w[4 * 4 * 11];
    // keyExpansion(key, w);

    for(int i = 0; i < 16; i++) {
        out[i] = in[i];
    }

    addRoundKey(out, w + 10 * 16);
    for(int i = 9; i > 0; i--) {
        invShiftRows(out);
        invSubBytes(out);
        addRoundKey(out, w + i * 16);
        invMixColumns(out);
    }
    invShiftRows(out);
    invSubBytes(out);
    addRoundKey(out, w);

    #ifdef DEBUG
    std::cout << endl << "---DecryptBlock Begin---" << endl;
    for(int i = 0; i < 16; i++) {
        std::cout << std::hex << (int)out[i] << " ";
    }
    std::cout << endl << "---DecryptBlock End---" << endl;
    #endif
}

void encryptECB(unsigned char* in, unsigned char* out, unsigned char* key, int size) {
    int padding = 16 - size % 16;
    unsigned char last_block[16];
    
    int i;
    for(i = 0; i < 16-padding; i++) {
        last_block[i] = in[size - 16 + padding + i];
    }
    for(; i < 16; i++) {
        last_block[i] = padding;
    }

    unsigned char w[4 * 4 * 11];
    keyExpansion(key, w);

    int p;
    int full_block_size = size - size % 16;
    for(p = 0; p < full_block_size; p += 16) {
        encryptBlock(in + p, out + p, w);
    }
    encryptBlock(last_block, out + p, w);

    #ifdef DEBUG
    std::cout << "Padding: " << padding << endl;
    for(i = 0; i < 16; i++) {
        std::cout << std::hex << (int)last_block[i] << " ";
    }
    std::cout << endl << "---EncryptECB Begin---" << endl;
    for(int i = 0; i < size + padding; i++) {
        std::cout << std::hex << (int)out[i] << " ";
    }
    std::cout << endl << "---EncryptECB End---" << endl;
    #endif
}

void decryptECB(unsigned char* in, unsigned char* out, unsigned char* key, int size) {
    unsigned char w[4 * 4 * 11];
    keyExpansion(key, w);

    int full_block_size = size - size % 16;
    int p;
    for(p = 0; p < full_block_size; p += 16) {
        decryptBlock(in + p, out + p, w);
    }
    decryptBlock(in + p, out + p, w);
    int padding = out[size];
    assert((padding + size) % 16 == 0);

    #ifdef DEBUG
    std::cout << endl << "---DecryptECB Begin---" << endl;
    std::cout << "Padding: " << std::dec << padding << endl;
    for(int i = 0; i < size + padding; i++) {
        std::cout << std::hex << (int)out[i] << " ";
    }
    std::cout << endl << "---DecryptECB End---" << endl;
    #endif
}

void encryptCBC(unsigned char* in, unsigned char* out, unsigned char* key, unsigned char* iv, int size) {
    int padding = 16 - size % 16;
    unsigned char last_block[16];

    int i;
    for(i = 0; i < 16-padding; i++) {
        last_block[i] = in[size - 16 + padding + i];
    }
    for(; i < 16; i++) {
        last_block[i] = padding;
    }

    unsigned char w[4 * 4 * 11];
    keyExpansion(key, w);

    int p = 0;
    int full_block_size = size - size % 16;

    if(p < full_block_size) {
        #pragma GCC unroll 16
        for(int j = 0; j < 16; j++) {
            in[j] ^= iv[j];
        }
        encryptBlock(in, out, w);
        p += 16;
    }
    for(; p < full_block_size; p += 16) {
        #pragma GCC unroll 16
        for(int j = 0; j < 16; j++) {
            in[p + j] ^= out[p + j - 16];
        }
        encryptBlock(in + p, out + p, w);
    }
    if(full_block_size) {
        #pragma GCC unroll 16
        for(int j = 0; j < 16; j++) {
            last_block[j] ^= out[p + j - 16];
        }
    }
    else {
        #pragma GCC unroll 16
        for(int j = 0; j < 16; j++) {
            last_block[j] ^= iv[j];
        }
    }
    encryptBlock(last_block, out + p, w);

    #ifdef DEBUG
    std::cout << "Padding: " << padding << endl;
    // for(i = 0; i < 16; i++) {
    //     std::cout << std::hex << (int)last_block[i] << " ";
    // }
    std::cout << endl << "---EncryptCBC Begin---" << endl;
    for(int i = 0; i < size + padding; i++) {
        std::cout << std::hex << (int)out[i] << " ";
    }
    std::cout << endl << "---EncryptCBC End---" << endl;
    #endif
}

void decryptCBC(unsigned char* in, unsigned char* out, unsigned char* key, unsigned char* iv, int size) {
    unsigned char w[4 * 4 * 11];
    keyExpansion(key, w);

    int full_block_size = size - size % 16;
    int p = 0;
    if(p < full_block_size) {
        decryptBlock(in, out, w);
        #pragma GCC unroll 16
        for(int j = 0; j < 16; j++) {
            out[j] ^= iv[j];
        }
        p += 16;
    }
    for(; p < full_block_size; p += 16) {
        decryptBlock(in + p, out + p, w);
        #pragma GCC unroll 16
        for(int j = 0; j < 16; j++) {
            out[p + j] ^= in[p - 16 + j];
        }
    }
    decryptBlock(in + p, out + p, w);
    if(full_block_size) {
        #pragma GCC unroll 16
        for(int j = 0; j < 16; j++) {
            out[p + j] ^= in[p - 16 + j];
        }
    }
    else {
        #pragma GCC unroll 16
        for(int j = 0; j < 16; j++) {
            out[j] ^= iv[j];    
        }
    }

    int padding = out[size];
    #ifdef DEBUG
    std::cout << "size: " << size << endl;
    std::cout << "Padding: " << std::dec << padding << endl;
    #endif
    assert((padding + size) % 16 == 0);

    #ifdef DEBUG
    std::cout << endl << "---DecryptCBC Begin---" << endl;
    for(int i = 0; i < size + padding; i++) {
        std::cout << std::hex << (int)out[i] << " ";
    }
    std::cout << endl << "---DecryptCBC End---" << endl;
    #endif
}

void encryptCTR(unsigned char* in, unsigned char* out, unsigned char* key, unsigned char* nonce, int size) {
    // CTR can be done without padding
    // nonce is 8-byte
    unsigned char w[4 * 4 * 11];
    keyExpansion(key, w);

    unsigned char counter[16] = {};
    for(int i = 0; i < 8; i++) {
        counter[i] = nonce[i];
    }

    #pragma omp parallel for firstprivate(counter)
    for(int p = 0; p < size; p += 16) {
        int idx = 15;
        int ctr = p >> 4;
        while(ctr) {
            counter[idx] = ctr & 0xff;
            ctr >>= 8;
            idx--;
        }
        // for(int i = 0; i < 16; i++) {
        //     std::cout << std::hex << (int) counter[i] << " ";
        // }
        // std::cout << endl;
        encryptBlock(counter, out + p, w);
        #pragma GCC unroll 16
        for(int i = 0; i < 16; i++) {
            out[p + i] ^= in[p + i];
        }
        // for(int i = 15; i >= 0; i--) {
        //     counter[i]++;
        //     if(counter[i]) {
        //         break;
        //     }
        // }
    }

    #ifdef DEBUG
    std::cout << endl << "---EncryptCTR Begin---" << endl;
    for(int i = 0; i < size; i++) {
        std::cout << std::hex << (int)out[i] << " ";
    }
    std::cout << endl << "---EncryptCTR End---" << endl;
    #endif
}
