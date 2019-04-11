//
// Created by haines on 3/27/19.
//

#include "Decompression.h"
#include "Coder.h"


Decompression::Decompression() {

}

Decompression::Decompression(string file_name, bool auto_load) {
    this->file_name = file_name;
    if (auto_load) load_data(file_name);
}


void Decompression::load_data(string file_name) {
    cout << "fuck!" << endl;
    ifstream fin;
    fin.open(file_name, ios::binary | ios::in);
    fin.read((char *)&raw_size, sizeof(raw_size));
    out_stream.load(&fin);
    fin.close();
    cout << "load bin data done!" << endl;
}

void Decompression::do_decompression() {
    adpac_decoder("result2.bin");
    cout << "arithmetic decoede done!" << endl;
    read_bitstream();
    cout << "read bitstream done!" << endl;

    dequantization();
    cout << "dequantization done!" << endl;

    svd_decode();
    cout << "scd_decode done!" << endl;

    delta_decode();
    cout << "delta docode done!" << endl;

    lcf_reconstruction();
    cout << "lcf reconstruction done!" << endl;

    reconstruction();
    cout << "reconstruction done!" << endl;

}

void Decompression::adpac_decoder(string file_name) {
    in_stream.resize_capacity(raw_size);
    entropy_coder.decode(raw_size, &out_stream, &in_stream);
    ifstream fin;
    int alphabet, result_size, data_size;

    fin.open(file_name, ios::binary | ios::in);
    fin.read(reinterpret_cast<char *>(&alphabet), sizeof(alphabet));
    fin.read(reinterpret_cast<char *>(&result_size), sizeof(result_size));
    fin.read(reinterpret_cast<char *>(&data_size), sizeof(data_size));

    unsigned char* result = (unsigned char*)malloc(result_size * sizeof(unsigned char));
    fin.read(reinterpret_cast<char *>(result), result_size);
    fin.close();

    //Decoding
    int test_size = data_size;
    int* test = (int*)malloc(test_size*sizeof(int));
    test_size = Decoding(result, alphabet, test);
    char *tt = (char *)malloc(test_size);
    for (uint32_t i = 0 ; i != test_size; ++i) tt[i] = (char)test[i];
    ////////////////////////////////////////////////////////////////////////////////
    cout << "test size: " << test_size << endl;
    bstream.load(tt);
    cout << "bstream: " << bstream.query_byte_occupancy() << endl;

//
//
//
//    //allocate memory and read data
//    int data_size = getFileSize(tmp_file.c_str());
//    if (data_size < 0) {
//        printf("File not found\n");
//        return;
//    }
//    unsigned char* data = (unsigned char*)malloc(data_size * sizeof(unsigned char));
//    FILE* f = fopen(tmp_file.c_str() , "rb");
//    fread(data, 1, data_size, f);
//    fclose(f);
//
//    int alphabet = 256;
//
//    ////////////////////////////////////////////////////////////////////////////////
//    //Encoding
//    int result_size = data_size*2 + 1000;
//    unsigned char* result = (unsigned char*)malloc(result_size * sizeof(unsigned char));
//    result_size = Encode(data, data_size, alphabet, result);
//    ////////////////////////////////////////////////////////////////////////////////
////    remove(tmp_file.c_str());
//
//    ofstream fout;
//    fout.open("result2.bin", ios::binary | ios::out | ios::trunc);
//    fout.write(reinterpret_cast<char *>(result), result_size);
//    fout.close();
//
//    cout << "out: " << out_stream.query_occupancy() << endl;
//    cout << "out2: " << result_size * 8 << endl;

}

void Decompression::read_bitstream() {
    // read bitstream head
    evx::uint32 head_size = sizeof(BitstreamHead);
    in_stream.read_bytes(&bithead, &head_size);
    cout << "read bithead, " << sizeof(BitstreamHead) <<" bytes." << endl;

    BitstreamHead bhead;
    head_size = sizeof(BitstreamHead);
    bstream.read_bytes(&bhead, &head_size);


    // read seeds
    for (uint32_t i = 0 ; i != bithead.seed_num * 9; ++i) {
        int32_t t = 0;
        in_stream.read_bits(&t, reinterpret_cast<evx::uint32 *>(&bithead.seed_bits));
        seed_triangle_points.push_back(static_cast<double_t >(t));
    }
    cout << "seeds: " << in_stream.query_occupancy() << endl;

    // read delta
    delta = Matrix(bithead.seed_num * 9, bithead.frame_num - 1);
    for (uint32_t i = 0 ; i != delta.Getrows(); ++i) {
        for (uint32_t j = 0 ; j != delta.Getcols(); ++j) {
            int32_t t = 0;
            in_stream.read_bits(&t, reinterpret_cast<evx::uint32 *>(&bithead.delta_bits));
            delta(i, j) = static_cast<double_t >(t);
        }
    }
    cout << "delta: " << in_stream.query_occupancy() << endl;

    // read U/C information
    uint32_t cluster_num = bithead.seed_num;
    UL.resize(cluster_num);
    CL.resize(cluster_num);
    Urow_array.resize(cluster_num);
    Ucol_array.resize(cluster_num);
    Crow_array.resize(cluster_num);
    Ccol_array.resize(cluster_num);
    evx::uint32 lsize = cluster_num * sizeof(double_t);
    evx::uint32 msize = cluster_num * sizeof(int);

    in_stream.read_bytes(UL.data(), &lsize);
    in_stream.read_bytes(CL.data(), &lsize);
    in_stream.read_bytes(Urow_array.data(), &msize);
    in_stream.read_bytes(Ucol_array.data(), &msize);
    in_stream.read_bytes(Crow_array.data(), &msize);
    in_stream.read_bytes(Ccol_array.data(), &msize);
    cout << "uc info: " << in_stream.query_occupancy() << endl;


    // read U
    U.resize(bithead.seed_num);
    for (uint32_t i = 0 ; i != U.size(); ++i) {
        U[i] = Matrix(Urow_array[i], Ucol_array[i]);
        for (uint32_t m = 0; m != U[i].Getrows(); ++m)
            for (uint32_t n = 0; n != U[i].Getcols(); ++n) {
                int32_t t = 0;
                in_stream.read_bits(&t, reinterpret_cast<evx::uint32 *>(&bithead.pca_bits));
                U[i](m, n) = static_cast<double>(t);
            }
    }
    cout << "U: " << in_stream.query_occupancy() << endl;


    // read C
    C.resize(bithead.seed_num);
    for (uint32_t i = 0 ; i != C.size(); ++i) {
        C[i] = Matrix(Crow_array[i], Ccol_array[i]);
        for (uint32_t m = 0; m != C[i].Getrows(); ++m)
            for (uint32_t n = 0; n != C[i].Getcols(); ++n) {
                int32_t t = 0;
                in_stream.read_bits(&t, reinterpret_cast<evx::uint32 *>(&bithead.pca_bits));
                C[i](m, n) = static_cast<double>(t);
            }
    }
    cout << "C: " << in_stream.query_occupancy() << endl;

    uint32_t size = sizeof(uint32_t);
    uint32_t t = 0;
    // read clusters points index
    cluster_points_index.resize(bithead.seed_num);
    for (uint32_t i = 0; i != bithead.seed_num; ++i) {
        for (uint32_t j = 0 ; j != Ucol_array[i] / 3 ; ++j) {
            in_stream.read_bytes(&t, &size);
            cluster_points_index[i].push_back(t);
        }
    }

    // read seed triangle point index
    for (uint32_t i = 0 ; i != bithead.seed_num * 3; ++i) {
        in_stream.read_bytes(&t, &size);
        seed_points_index.push_back(t);
    }
}

void Decompression::dequantization() {
    for (uint32_t i = 0 ; i != seed_triangle_points.size() ; ++i)
        seed_triangle_points[i] = (seed_triangle_points[i] * bithead.seed_length) / pow(2, bithead.seed_bits);
    delta.dequantization(bithead.delta_bits, bithead.delta_length);
    for (uint32_t i = 0 ; i != U.size(); ++i) {
        U[i].dequantization(bithead.pca_bits, UL[i]);
        C[i].dequantization(bithead.pca_bits, CL[i]);
    }
}

void Decompression::svd_decode() {
    for (uint32_t i = 0; i != U.size(); ++i) {
        data.push_back(U[i].trans() * C[i]);
    }
}

void Decompression::delta_decode() {
    seed_data = Matrix(bithead.seed_num * 9, bithead.frame_num);
    for (uint32_t i = 0; i != bithead.frame_num; ++i) {
        for (uint32_t j = 0; j != seed_data.Getrows(); ++j) {
            if (i == 0)
                seed_data(j, i) = seed_triangle_points[j];
            else{
               seed_data(j, i) = seed_data(j, i - 1) + delta(j, i - 1);
            }
        }
    }
}

void Decompression::reconstruction() {
    num_vertices = seed_points_index.size();
    for (uint32_t i = 0 ; i != cluster_points_index.size(); ++i)
        num_vertices += cluster_points_index[i].size();
    cout << "vertices number: " << num_vertices << endl;
    vertices.resize(num_vertices);

    Point p, pw;
    for (uint32_t i = 0; i != cluster_points_index.size(); ++i) {
        auto& cluster = cluster_points_index[i];
        for (uint32_t j = 0 ; j != cluster.size(); ++j) {
            for (uint32_t f = 0 ; f != bithead.frame_num; ++f) {
                LC lc = lcfs[f].lcs[i];
                p.x = data[i](3 * j, f);
                p.y = data[i](3 * j + 1, f);
                p.z = data[i](3 * j + 2, f);

                pw = lc.local2world(p);
                vertices[cluster[j]].push_back(pw);
            }
        }
    }

    for (uint32_t i = 0; i != seed_points_index.size(); ++i) {
        for (uint32_t j = 0 ; j != bithead.frame_num; ++j) {
            p.x = seed_data(3 * i, j);
            p.y = seed_data(3 * i + 1, j);
            p.z = seed_data(3 * i + 2, j);
            vertices[seed_points_index[i]].push_back(p);
        }
    }

}

void Decompression::lcf_reconstruction(int frame_idx) {
    LCF lcf;
    Point p1, p2, p3;
    lcf.frame_idx = frame_idx;
    for (uint32_t i = 0; i != seed_data.Getrows(); i += 9) {
        p1.x = seed_data(i, frame_idx);
        p1.y = seed_data(i + 1, frame_idx);
        p1.z = seed_data(i + 2, frame_idx);
        p2.x = seed_data(i + 3, frame_idx);
        p2.y = seed_data(i + 4, frame_idx);
        p2.z = seed_data(i + 5, frame_idx);
        p3.x = seed_data(i + 6, frame_idx);
        p3.y = seed_data(i + 7, frame_idx);
        p3.z = seed_data(i + 8, frame_idx);

        LC lc;

        lc.o = (p1 + p2) / 2;


        Vector3 p12 = p2 - p1;
        Vector3 p13 = p3 - p1;
        Vector3 normal = p12 ^ p13;
        Vector3 vy = p12 ^ normal;

        lc.ix = p12 / vec_length(p12);
        lc.iy = vy / vec_length(vy);
        lc.iz = normal / vec_length(normal);

        Matrix33 m;

        m << lc.ix.x, lc.iy.x, lc.iz.x,
                lc.ix.y, lc.iy.y, lc.iz.y,
                lc.ix.z, lc.iy.z, lc.iz.z;

        lc.affline_matrix = m.inverse();
        lc.restore_matrix = m;

        lcf.lcs.push_back(lc);
    }
    lcfs.push_back(lcf);
}

void Decompression::lcf_reconstruction() {
    for (uint32_t i = 0; i != bithead.frame_num; ++i)
        lcf_reconstruction(i);
}

void Decompression::save_objs(string data_dir) {

}



