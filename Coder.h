//
// Created by haines on 4/8/19.
//

#ifndef ANIMATION_COMPRESSION_CODER_H
#define ANIMATION_COMPRESSION_CODER_H

// RunCoder15.cpp Adaptive order 1.5 range coder
// (C) 2010, Andrew Polar under GPL ver. 3.
// Released  Feb, 2010.
//
//   LICENSE
//
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License as
//   published by the Free Software Foundation; either version 3 of
//   the License, or (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful, but
//   WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   General Public License for more details at
//   Visit <http://www.gnu.org/copyleft/gpl.html>.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

//By commenting out this flag you switch from order 1.5 range coder to order 0.
#define USE_PREDICTOR

//Globals
//The following two constants can be changed by someone who understands the basic
//algorithm. When alphabet is large MAX_CUMULATIVE_FREQUENCY_BIT_SIZE should be larger.
#define MANTISSA_BIT_SIZE 16                 //defines computational precision
#define MAX_CUMULATIVE_FREQUENCY_BIT_SIZE 14 //defines precision for expressing frequencies
#define LENGTH_OF_ADAPTATION_BUFFER 4096     //defines number of symbols after which we recalculate frequencies

extern const int			MAX_INDEX ;
extern const unsigned char Length_of_context_buffer; //when changed function get_selected_context should be changed
extern const int			context_size;

extern unsigned			current_byte;
extern unsigned long long  operation_buffer64;
extern unsigned long long  bit_holding_buffer64;
extern unsigned char       bit_counter;
extern unsigned int        operation_buffer32;
extern int                 MANTISSA, EXPONENT;
extern unsigned char*      global_pointer_to_data_buffer;
extern long long			overflow_indicator;

static __inline int get_selected_context(int memory_buffer) {
    return (memory_buffer & 0xFF) | ((memory_buffer & 0xF000)>>4);
}

static __inline unsigned char GetBitLength(int N) {
    unsigned char len = 0;
    while ((N >>= 1) > 0)
        ++len;
    return len + 1;
}

static __inline void SafeProduct(int F, unsigned char L) {
    MANTISSA *= F;
    if ((MANTISSA & (1 << (MANTISSA_BIT_SIZE + L - 1))) > 0) {
        EXPONENT = L;
    }
    else {
        EXPONENT = L - 1;
    }
    MANTISSA >>= EXPONENT;
}

static __inline void readBits(int bits) {
    if (bits == 0) {
        return;
    }

    if ((bit_counter - bits) < 0) {
        bit_holding_buffer64 <<= 8;
        bit_holding_buffer64 |= global_pointer_to_data_buffer[current_byte+0];
        bit_holding_buffer64 <<= 8;
        bit_holding_buffer64 |= global_pointer_to_data_buffer[current_byte+1];
        bit_holding_buffer64 <<= 8;
        bit_holding_buffer64 |= global_pointer_to_data_buffer[current_byte+2];
        bit_holding_buffer64 <<= 8;
        bit_holding_buffer64 |= global_pointer_to_data_buffer[current_byte+3];
        current_byte += 4;
        bit_counter += 32;
    }
    operation_buffer64 <<= bits;
    operation_buffer64 |= (bit_holding_buffer64 >> (bit_counter - bits)) & (0xffffffff >> (32 - bits));
    bit_counter -= bits;
}

static __inline void writeBits(int bits) {
    if (bits == 0) {
        return;
    }

    bit_holding_buffer64 <<= bits;
    bit_holding_buffer64 |= (operation_buffer64 >> (64-bits));
    bit_counter += bits;
    if (bit_counter > 31) {
        unsigned char offset = bit_counter - 32;
        global_pointer_to_data_buffer[current_byte+0] = (unsigned char)(bit_holding_buffer64 >> (offset + 24));
        global_pointer_to_data_buffer[current_byte+1] = (unsigned char)(bit_holding_buffer64 >> (offset + 16));
        global_pointer_to_data_buffer[current_byte+2] = (unsigned char)(bit_holding_buffer64 >> (offset + 8));
        global_pointer_to_data_buffer[current_byte+3] = (unsigned char)(bit_holding_buffer64 >> offset);
        bit_counter -= 32;
        current_byte += 4;
    }
}

int getFileSize(const char* fileName);

#ifdef USE_PREDICTOR

//adaptive order 1.5 predictor
class Predictor {
public:
    Predictor(int History, int Prediction);
    ~Predictor();
    void SetUpInitial();
    void UpdateTables(int address, unsigned char value);
    unsigned char GetSymbol(int address, unsigned char value);
    unsigned char InvertSymbol(int address, unsigned char value);
private:
    unsigned char** frequency;
    unsigned char**  symbol;
    unsigned char**  inverse;
    int history_size, predict_size;
};
#endif //USE_PREDICTOR


//Adaptive oder 0 range coder
class RunCoder {
public:
    RunCoder();
    ~RunCoder();
    bool Initialize(int size);
    void Renormalize(bool isDecoding);
    void ProcessSymbol(int value);
    void FlushRemainedBuffer();
    int  DecodeSymbol();
    void UpdateRawFrequency(int value, bool isDecoding);
    void SetDefaultFrequencies();
private:
    int* rawfrequency;
    int* normalfrequency;
    int* cumulative;
    int* symbol;
    unsigned char* bit_freq_len;
    void FreeResources();
    int  alphabetsize;
    int  memorylimit;
    int nCounter;
};

template <class T>
int Encode(T* data, int data_size, int alphabet, unsigned char* result) {
    clock_t start_encoding = clock();
    printf("Start encoding\n");

    global_pointer_to_data_buffer = result;
    RunCoder* encoder = new RunCoder();
    if (encoder->Initialize(alphabet) == false) {
        printf("Initialization failed\n");
        exit(1);
    }
    encoder->SetDefaultFrequencies();
    encoder->Renormalize(false);

#ifdef USE_PREDICTOR
    Predictor* pp = new Predictor(context_size, alphabet);
#endif

    int context = 0;
    unsigned char symbol;
    for (int i=0; i<data_size; ++i) {

#ifdef USE_PREDICTOR
        //alphabet replacement
        symbol = pp->GetSymbol(context, data[i]);
        pp->UpdateTables(context, data[i]);
        context <<= 8;
        context |= data[i];
        context = get_selected_context(context);
        //end
#endif

        //here we check for buffer overflow
        while (((operation_buffer64 << (MAX_CUMULATIVE_FREQUENCY_BIT_SIZE - EXPONENT)) & overflow_indicator) == overflow_indicator) {
            //here we correct buffer overflow
            printf("Buffer overflow is corrected at %d\n", i);
            encoder->ProcessSymbol(-2);
        }

#ifdef USE_PREDICTOR
        encoder->ProcessSymbol(symbol);
        encoder->UpdateRawFrequency(symbol, false);
#else
        encoder->ProcessSymbol(data[i]);
		encoder->UpdateRawFrequency(data[i], false);
#endif
    }
    encoder->FlushRemainedBuffer();
#ifdef USE_PREDICTOR
    delete pp;
#endif
    delete encoder;

    clock_t end_encoding = clock();
    printf("End encoding, time %2.3f s., size   %d \n", (double)(end_encoding - start_encoding)/CLOCKS_PER_SEC, current_byte);
    return current_byte;
}

template <class T>
int Decoding(unsigned char* code, int alphabet, T* decode) {
    clock_t start_decoding = clock();
    printf("Start decoding\n");
    global_pointer_to_data_buffer = code;

    RunCoder* decoder = new RunCoder();
    if (decoder->Initialize(alphabet) == false) {
        printf("Unexpected bug\n");
        exit(1);
    }
    decoder->SetDefaultFrequencies();
    decoder->Renormalize(true);

#ifdef USE_PREDICTOR
    Predictor* pp = new Predictor(context_size, 256);
#endif

    int cnt = 0;
    int value;
    int context = 0;
    unsigned char symbol;
    bool isOK = false;
    readBits(32);
    readBits(32);
    while (true) {
        value = decoder->DecodeSymbol();
        if (value == 1) {
            isOK = true;
            break;
        }
        if (value == 0) {
            printf("Zero flag recognized at %d\n", cnt);
        }
        else {
            value -=2;
            decoder->UpdateRawFrequency(value, true);
#ifdef USE_PREDICTOR
            //revert alphabet
            symbol = pp->InvertSymbol(context, value);
            pp->UpdateTables(context, symbol);
            //end
            context <<= 8;
            context |= symbol;
            context = get_selected_context(context);
            decode[cnt++] = symbol;
#else
            decode[cnt++] = value;
#endif
        }
    }
#ifdef USE_PREDICTOR
    delete pp;
#endif
    delete decoder;

    clock_t end_decoding = clock();
    printf("End decoding, time %2.3f s.\n", (double)(end_decoding - start_decoding)/CLOCKS_PER_SEC);
    if (isOK == true)
        printf("Termination flag is recognized\n");
    return cnt;
}


#endif //ANIMATION_COMPRESSION_CODER_H
