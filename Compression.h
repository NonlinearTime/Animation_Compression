//
// Created by haines on 3/27/19.
//

#ifndef ANIMATION_COMPRESSION_COMPRESSION_H
#define ANIMATION_COMPRESSION_COMPRESSION_H

#include <iostream>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include "Cluster.h"
#include "WavefrontObjFile.h"
#include "Eigen.h"
#include "Matrix.h"
#include "bitstream.h"
#include "cabac.h"



using namespace std;

using DeltaRow = vector<Point >;


struct BitstreamHead {
    int seed_num;
    int frame_num;
    int face_num;
    int seed_bits;
    double seed_length;
    double seed_low;
    int delta_bits;
    double delta_length;
    double delta_low;
    int pca_bits;
    uint32_t raw_size;
};


class Compression {
    Cluster clst;
//    vector<DeltaRow > delta;
    Matrix delta;
    vector<double > seed_triangle_points;
    vector<Matrix > U, C;
    vector<double > UL, CL;
    vector<double > Ul, Cl;
    vector<int > Urow_array, Crow_array, Ucol_array, Ccol_array;
    BitstreamHead bithead;
    evx::bitstream in_stream, out_stream;
    evx::entropy_coder entropy_coder;
    vector<uint32_t > seed_points_index;
    vector<cluster > cluster_points_index;
    uint32_t raw_size;

public:
    Compression () {}

    ~Compression() {}

    void delta_encode();

    void pca_encode();
    Matrix pca_decode(Matrix C, Matrix V);
    void svd_encode();
    Matrix svd_decode(Matrix C, Matrix U);

    void set_cluster_num(uint32_t N);

    void do_compression();
    void do_compression(uint32_t cluster_num);

    void save_data(string file_name);

    bool load_data(string data_dir, int first_frame, int num_frames);

    void rate_distortion();

    void rd_optimization();

    void quantitzation(int pca_bits, int seed_bits, int delta_bits);

    void write_bitstream();

    double mse_error(Matrix p, Matrix q);

    void adpac_encoder();

private:

    void pca_per_cluster(Matrix m, int N);
    void svd_per_cluster(Matrix m, int N);

};


#endif //ANIMATION_COMPRESSION_COMPRESSION_H
