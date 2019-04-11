#include <iostream>
#include <fstream>
#include "WavefrontObjFile.h"
#include "Cluster.h"
#include "Compression.h"
#include "Decompression.h"
#include "bitstream.h"
#include "Coder.h"

int main() {
    std::cout << "Hello, World!" << std::endl;
    Point3D<int> X(1,2, 3.1);
    Point3D<int> Y(2,2,2);
    X += Y;
    X.Print();
    cout << X;
    int x = 1;
    int &y = x;
    y = 2;
    cout << x << " " << y << endl;
//    char t[20];
//    int a , b;
//    scanf("%s ", t);
//    cout << t << endl;
//    cout << a << " " << b << endl;
    WavefrontObjFile<float_t > wObj;
    char fmt[256]  = "/home/haines/Documents/compression/Reconstructed_mesh/%06d.obj";
    char buf[256];
    sprintf(buf, fmt, 20);
    cout << buf << endl;
//    for (int i = 750; i <= 870; ++i) {
//        sprintf(buf, fmt, i);
//        wObj.load_obj(buf);
//    }
    Matrix33 m;
    m << 1,2,3,
         4,5,6,
         7,8,9;
    cout << m << endl;

    vector<cluster > clusters;
    cout << clusters.size() << endl;
    clusters.resize(10);
    cout << clusters[1].size() << endl;

    Obj<float_t > obj, next_obj;
    WavefrontObjFile<float_t > tObj;
    tObj.load_obj("/home/haines/Documents/compression/Reconstructed_mesh/000000.obj");
    obj = move(tObj.get_obj());

    fstream f_obj;
    f_obj.open("obj_poses.txt", ios::out | ios::trunc);
    f_obj << 44 << " " << obj.vertices.size() << endl;

    string tmp = "hello";
    string w = "world!\n";
    tmp.append(w);
    cout << tmp ;

//    for (uint32_t i = 1; i != 45; ++i) {
//        sprintf(buf, fmt, i);
//        wObj.load_obj(buf);
//        next_obj = wObj.get_obj();
//        for (uint32_t j = 0 ; j != obj.vertices.size() ; ++j) {
//            f_obj << next_obj.vertices[j] - obj.vertices[j] << " ";
//        }
//        f_obj << endl;
//        obj = wObj.get_obj();
//    }
//
//    f_obj.close();

//
//    ofstream fout;
//
//    char c = 'c';
//    fout.open("test.txt", ios::trunc | ios::out | ios::binary);
//    fout.write(&c, 1);
//    fout.close();
//
//    Compression compression;
//    compression.load_data("/home/haines/Documents/compression/Reconstructed_mesh", 0, 40);
//    compression.set_cluster_num(200);
//    compression.do_compression();
//    compression.save_data("result.bin");

//    BitstreamHead bhead;
//    bhead.seed_num = 1;
//    bhead.frame_num = 2;
//    bhead.seed_bits = 3;
//    bhead.seed_length = 0.9;
//    bhead.delta_bits = 4;
//    bhead.delta_length = 1.8;
//    bhead.pca_bits = 5;
//    bhead.raw_size = 6;
//
//    int data_size = 48;
//    int alphabet = 256;
//    int result_size = data_size*2 + 1000;
//    unsigned char* result = (unsigned char*)malloc(result_size * sizeof(unsigned char));
//    result_size = Encode((char *)&bhead, data_size, alphabet, result);
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    //Decoding
//    int test_size = data_size;
//    int* test = (int*)malloc(test_size*sizeof(int));
//    test_size = Decoding(result, alphabet, test);
//    char *tt = (char *)malloc(test_size);
//    for (uint32_t i = 0 ; i != test_size; ++i) tt[i] = (char)test[i];
//
//    BitstreamHead t;
//    memcpy(&t, tt, data_size);
//
//    evx::bitstream in(1024), out(1024), test(1024);
//    in.write_bytes(&bhead, 48);


//    BitstreamHead bbhead;
//    uint32_t size = 48;
////    in.read_bytes(&bbhead, &size);
//    uint32_t raw_size = in.query_occupancy();
//
//    evx::entropy_coder coder;
//    coder.encode(&in, &out);
//    coder.decode(raw_size, &out, &test);
//
//    BitstreamHead thead;
//    size = 48;
//    test.read_bytes(&thead, &size);


    Decompression decompression;
    decompression.load_data("result.bin");
    decompression.do_decompression();

    return 0;
}