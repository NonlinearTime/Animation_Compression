//
// Created by haines on 4/10/19.
//

#include "Coder.h"
//
const int			MAX_INDEX = (1 << MAX_CUMULATIVE_FREQUENCY_BIT_SIZE)-1;
const unsigned char Length_of_context_buffer = 12; //when changed function get_selected_context should be changed
const int			context_size = 1<<Length_of_context_buffer;

unsigned			current_byte = 0;
unsigned long long  operation_buffer64 = 0;
unsigned long long  bit_holding_buffer64 = 0;
unsigned char       bit_counter = 0;
unsigned int        operation_buffer32 = 0;
int                 MANTISSA = 0, EXPONENT = 0;
unsigned char*      global_pointer_to_data_buffer = NULL;
long long			overflow_indicator = 0xffffffffffffffff << (MAX_CUMULATIVE_FREQUENCY_BIT_SIZE + MANTISSA_BIT_SIZE);

#ifdef USE_PREDICTOR

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

int getFileSize(const char* fileName) {

    FILE* f;
    if ((f = fopen(fileName, "rb")) == NULL)
        return -1;

    fseek(f, 0, SEEK_END);
    int bytes = ftell(f);
    fclose(f);

    return bytes;
}