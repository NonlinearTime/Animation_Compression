#include <iostream>
#include <fstream>
#include "WavefrontObjFile.h"
#include "Cluster.h"

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


    ofstream fout;

    char c = 'c';
    fout.open("test.txt", ios::trunc | ios::out | ios::binary);
    fout.write(&c, 1);
    fout.close();


    return 0;
}