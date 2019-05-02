//
// Created by haines on 3/27/19.
//

#ifndef ANIMATION_COMPRESSION_DECOMPRESSION_H
#define ANIMATION_COMPRESSION_DECOMPRESSION_H

#include <cstring>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <map>
#include <algorithm>
#include <cstdlib>
#include "Compression.h"

using namespace std;



class Decompression {
//    vector<DeltaRow > delta;
    Matrix delta;
    vector<double > seed_triangle_points;
    vector<Matrix > U, C;
    vector<double > UL, CL;
    vector<double > Ul, Cl;
    vector<int > Urow_array, Crow_array, Ucol_array, Ccol_array;
    BitstreamHead bithead;
    evx::bitstream in_stream, out_stream, bstream;
    evx::entropy_coder entropy_coder;

    vector<uint32_t > seed_points_index;
    vector<cluster > cluster_points_index;

    vector<Matrix > data;
    Matrix seed_data;
    vector<vector<Point >> vertices;
    uint32_t num_vertices;

    vector<TriangleFace> faces;

    Frames frames;

    vector<LCF > lcfs;

    string file_name;

    uint32_t raw_size;

public:
    Decompression();
    Decompression(string file_name, bool auto_load = true);

    void do_decompression();
    void load_data(string file_name);
    void save_objs(string data_dir);

private:
    void adpac_decoder(string file_name);
    void read_bitstream();
    void dequantization();
    void svd_decode();
    void delta_decode();
    void lcf_reconstruction();
    void lcf_reconstruction(int frame_idx);
    void reconstruction();
};


#endif //ANIMATION_COMPRESSION_DECOMPRESSION_H
