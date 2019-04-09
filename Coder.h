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

const int			MAX_INDEX = (1 << MAX_CUMULATIVE_FREQUENCY_BIT_SIZE)-1;
const unsigned char Length_of_context_buffer = 12; //when changed function get_selected_context should be changed
const int			context_size = 1<<Length_of_context_buffer;

unsigned			current_byte;
unsigned long long  operation_buffer64;
unsigned long long  bit_holding_buffer64;
unsigned char       bit_counter;
unsigned int        operation_buffer32;
int                 MANTISSA, EXPONENT;
unsigned char*      global_pointer_to_data_buffer;
long long			overflow_indicator = 0xffffffffffffffff << (MAX_CUMULATIVE_FREQUENCY_BIT_SIZE + MANTISSA_BIT_SIZE);

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

Predictor::Predictor(int History, int Prediction) {
    history_size = History;
    predict_size = Prediction;

    frequency = (unsigned char**) malloc(history_size*sizeof(unsigned char*));
    symbol    = (unsigned char**) malloc(history_size*sizeof(unsigned char*));
    inverse   = (unsigned char**) malloc(history_size*sizeof(unsigned char*));

    for (int i=0; i<history_size; ++i) {
        frequency[i] = (unsigned char*) malloc(predict_size*sizeof(unsigned char));
        symbol[i]    = (unsigned char*) malloc(predict_size*sizeof(unsigned char));
        inverse[i]   = (unsigned char*) malloc(predict_size*sizeof(unsigned char));

        for (int j=0; j<predict_size; ++j) {
            frequency[i][j] = 0; symbol[i][j] = j; inverse[i][j] = j;
        }
    }
}

Predictor::~Predictor() {
    if (frequency) {
        for (int i=0; i<history_size; ++i) {
            if (frequency[i]) {
                free(frequency[i]);
                free(symbol[i]);
                free(inverse[i]);
            }
        }
        free(frequency);
        free(symbol);
        free(inverse);
    }
}

void Predictor::UpdateTables(int address, unsigned char value) {
    ++frequency[address][value];
    if (frequency[address][value] == 0xff) {
        for (int j=0; j<predict_size; ++j) {
            frequency[address][j] >>= 1;
        }
    }

    if (symbol[address][value] == 0) {
        return;
    }

    unsigned char  pos = symbol[address][value];
    unsigned int Freq = frequency[address][value];

    int npos = -1;
    for (int i=pos-1; i>=0; --i) {
        if (Freq > frequency[address][inverse[address][i]]) {
            npos = i;
        }
        else {
            break;
        }
    }
    if (npos < 0) {
        return;
    }

    unsigned char buf = symbol[address][inverse[address][pos]];
    symbol[address][inverse[address][pos]] = symbol[address][inverse[address][npos]];
    symbol[address][inverse[address][npos]] = buf;

    buf = inverse[address][pos];
    inverse[address][pos] = inverse[address][npos];
    inverse[address][npos] = buf;
}

unsigned char Predictor::GetSymbol(int address, unsigned char value) {
    return symbol[address][value];
}

unsigned char Predictor::InvertSymbol(int address, unsigned char value) {
    return inverse[address][value];
}
//End predictor
/////////////////////////////////////////////////////////////////////
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

RunCoder::RunCoder() {
    rawfrequency = 0;
    normalfrequency  = 0;
    cumulative = 0;
    symbol     = 0;
    bit_freq_len = 0;
}

RunCoder::~RunCoder() {
    FreeResources();
}

bool RunCoder::Initialize(int size) {
    if (rawfrequency == 0) {
        rawfrequency = (int*)malloc(size * sizeof(int));
        if (rawfrequency == 0)
            return false;
        memset(rawfrequency, 0x00, size * sizeof(int));
    }
    alphabetsize = size + 2;

    if (normalfrequency == 0) {
        normalfrequency = (int*)malloc(alphabetsize * sizeof(int));
        if (normalfrequency == 0)
            return false;
    }

    if (bit_freq_len == 0) {
        bit_freq_len = (unsigned char*)malloc(alphabetsize * sizeof(unsigned char));
        if (bit_freq_len == 0)
            return false;
    }

    if (cumulative == 0) {
        cumulative = (int*)malloc((alphabetsize+1)*sizeof(int)); //must be size+1
        if (cumulative == 0)
            return false;
    }

    if (symbol == 0) {
        symbol = (int*)malloc((1 << MAX_CUMULATIVE_FREQUENCY_BIT_SIZE) * sizeof(int));
        if (symbol == 0)
            return false;
    }

    operation_buffer64 = 0;
    bit_holding_buffer64 = 0;
    EXPONENT = 0;
    current_byte = 0;
    bit_counter = 0;
    operation_buffer32 = 0;
    MANTISSA = (1<<MANTISSA_BIT_SIZE) - 1;
    memorylimit = 0xfff;
    nCounter = 0;

    return true;
}

void RunCoder::FreeResources() {
    if (rawfrequency) {
        free(rawfrequency);
        rawfrequency = 0;
    }
    if (normalfrequency) {
        free(normalfrequency);
        normalfrequency = 0;
    }
    if (bit_freq_len) {
        free(bit_freq_len);
        bit_freq_len = 0;
    }
    if (cumulative) {
        free(cumulative);
        cumulative = 0;
    }
    if (symbol) {
        free(symbol);
        symbol = 0;
    }
}

void RunCoder::Renormalize(bool isDecoding) {
    normalfrequency[0] = 1; //symbol 0 is reserved for outputting bits and avoiding buffer overflow
    normalfrequency[1] = 1; //symbol 1 is reserved as end of data marker, when decoder reads 1 it stops processing.

    //other frequencies are copied to local buffer
    for (int i=2; i<alphabetsize; ++i) {
        normalfrequency[i] = rawfrequency[i-2];
    }

    //normalization of frequency
    int iTotal = 0;
    for (int i=0; i<alphabetsize; ++i)
        iTotal += normalfrequency[i];

    for (int i=2; i<alphabetsize; ++i) { //we start from 2 because 0 and 1 are reserved flags
        normalfrequency[i] = (normalfrequency[i] << MAX_CUMULATIVE_FREQUENCY_BIT_SIZE) / iTotal;
        if (normalfrequency[i] == 0)
            normalfrequency[i] = 1;
    }

    //recalculate Total
    iTotal = 0;
    for (int i=0; i<alphabetsize; ++i) {
        iTotal += normalfrequency[i];
    }

    if (iTotal > MAX_INDEX) {  //we need correction, iTotal should not be > MAX_INDEX.
        int i = 2;
        while (iTotal > MAX_INDEX) {
            if (normalfrequency[i] > 1) {
                --normalfrequency[i];
                --iTotal;
            }
            ++i;
            if (i == alphabetsize)
                i = 2;
        }
    }
    //end of normalization

    //we calculate bit lengthes of every frequency
    for (int i=0; i<alphabetsize; ++i) {
        bit_freq_len[i] = GetBitLength(normalfrequency[i]);
    }

    //making cumulative frequencies
    cumulative[0] = 0;
    for (int i=1; i<alphabetsize; ++i) {
        cumulative[i] = cumulative[i-1] + normalfrequency[i-1];
    }
    cumulative[alphabetsize] = (1 << MAX_CUMULATIVE_FREQUENCY_BIT_SIZE); //this is done on purpose

    if (isDecoding == false)
        return;

    //making symbol lookup table
    for (int k=0; k<alphabetsize; ++k) {
        for (int i=cumulative[k]; i<cumulative[k+1]; ++i) {
            symbol[i] = k;
        }
    }
}

void RunCoder::UpdateRawFrequency(int value, bool isDecoding) {
    ++rawfrequency[value];
    ++nCounter;
    if (rawfrequency[value] >= memorylimit) { //this is memory mechanism
        for (int i=0; i<alphabetsize-2; ++i) {
            rawfrequency[i] >>= 1;
        }
    }
    if (nCounter == LENGTH_OF_ADAPTATION_BUFFER) {
        Renormalize(isDecoding);
        nCounter = 0;
    }
}

void RunCoder::ProcessSymbol(int value) {
    writeBits(MAX_CUMULATIVE_FREQUENCY_BIT_SIZE - EXPONENT);
    operation_buffer64 <<= MAX_CUMULATIVE_FREQUENCY_BIT_SIZE - EXPONENT;
    operation_buffer64 += cumulative[value+2] * MANTISSA;
    SafeProduct(normalfrequency[value+2], bit_freq_len[value+2]);
}

void RunCoder::FlushRemainedBuffer() {
    ProcessSymbol(-1);

    //flushing 64 bits of operation_buffer64
    writeBits(32);
    operation_buffer64 <<= 32;
    writeBits(32);
    operation_buffer64 <<= 32;

    //we write some remained bits delayed by writing
    if (bit_counter > 0) {
        writeBits(32-bit_counter);
    }
}

int RunCoder::DecodeSymbol() {
    readBits(MAX_CUMULATIVE_FREQUENCY_BIT_SIZE - EXPONENT);
    int	ID = int(operation_buffer64/MANTISSA);
    if (ID > MAX_INDEX || ID < 0) {
        printf("Error in decoding, process terminated.\n");
        return 1;
    }
    operation_buffer64 -= MANTISSA * cumulative[symbol[ID]];
    SafeProduct(normalfrequency[symbol[ID]], bit_freq_len[symbol[ID]]);
    return symbol[ID];
}

void RunCoder::SetDefaultFrequencies() {
    for (int i=0; i<alphabetsize-2; ++i) {
        rawfrequency[i] = 1;
    }
#ifdef USE_PREDICTOR
    //this works with exponential distribution only
    rawfrequency[0] = 240; rawfrequency[1] = 220; rawfrequency[2] = 180; rawfrequency[3] = 120;
    rawfrequency[4] = 90;  rawfrequency[5] = 80;  rawfrequency[6] = 60;
#endif
}
//End RunCoder
/////////////////////////////////////////////////////

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

int getFileSize(const char* fileName) {

    FILE* f;
    if ((f = fopen(fileName, "rb")) == NULL)
        return -1;

    fseek(f, 0, SEEK_END);
    int bytes = ftell(f);
    fclose(f);

    return bytes;
}

#endif //ANIMATION_COMPRESSION_CODER_H
