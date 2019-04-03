//
// Created by haines on 3/27/19.
//

#include "Compression.h"

void Compression::save_data(char *file_name) {
    // write out_bitstream into bin file

    out_stream.save(file_name);
}

bool Compression::load_data(string data_dir, int first_frame, int num_frames) {
    if (data_dir[data_dir.length()] != '/') data_dir.append("/");
    string file_name = data_dir;
    char buf[256];
    file_name.append("%06d.obj");
    Frames frames;
    WavefrontObjFile<float_t > wObj;
    for (uint32_t i = 0 ; i != num_frames ; ++i) {
        sprintf(buf, file_name.c_str(), i);
        wObj.load_obj(buf);
        frames.push_back(wObj.get_obj());
    }
    clst.set_frames(frames);
    return true;
}

void Compression::do_compression() {
    clst.segmentation();
    delta_encode();
    svd_encode();
    quantitzation(16, 16, 16);
    write_bitstream();
    adpac_encoder();
}

void Compression::set_cluster_num(uint32_t N) {
    clst.set_clusters_num(N);
}

void Compression::do_compression(uint32_t cluster_num) {
    clst.set_clusters_num(cluster_num);
    do_compression();
}

void Compression::delta_encode() {
    auto frames = clst.get_frames();
    auto seed_triangles = clst.get_seed_triangles();

    Point p1, p2, p3;
    Point t1, t2, t3;
    delta = Matrix(9 * clst.get_clusters_num(), frames.size() - 1 );

    for (uint32_t i = 0 ; i != seed_triangles.size() ; ++i) {
        for (uint32_t j = 0 ; j != frames.size() ; ++j) {
            p1 = frames[j].vertices[seed_triangles[i].v[0]];
            p2 = frames[j].vertices[seed_triangles[i].v[1]];
            p3 = frames[j].vertices[seed_triangles[i].v[2]];
            if (j == 0) {
                seed_triangle_points.push_back(p1.x);
                seed_triangle_points.push_back(p1.y);
                seed_triangle_points.push_back(p1.z);
                seed_triangle_points.push_back(p2.x);
                seed_triangle_points.push_back(p2.y);
                seed_triangle_points.push_back(p2.z);
                seed_triangle_points.push_back(p3.x);
                seed_triangle_points.push_back(p3.y);
                seed_triangle_points.push_back(p3.z);
            } else {
                t1 = p1 - t1;
                t2 = p2 - t2;
                t3 = p3 - t3;
                delta(9 * i, j) = t1.x;
                delta(9 * i + 1, j) = t1.y;
                delta(9 * i + 2, j) = t1.z;
                delta(9 * i + 3, j) = t2.x;
                delta(9 * i + 4, j) = t2.y;
                delta(9 * i + 5, j) = t2.z;
                delta(9 * i + 6, j) = t3.x;
                delta(9 * i + 7, j) = t3.y;
                delta(9 * i + 8, j) = t3.z;
            }
            t1 = p1;
            t2 = p2;
            t3 = p3;
        }
    }

//    for (uint32_t i = 0 ; i != seed_triangles.size() ; ++i) {
//        DeltaRow dr;
//        for (uint32_t j = 0 ; j != frames.size() ; ++j) {
//            p1 = frames[j].vertices[seed_triangles[i].v[0]];
//            p2 = frames[j].vertices[seed_triangles[i].v[1]];
//            p3 = frames[j].vertices[seed_triangles[i].v[2]];
//            if (j == 0) {
//                dr.push_back(p1);
//                dr.push_back(p2);
//                dr.push_back(p3);
//            } else {
//                dr.push_back(p1 - t1);
//                dr.push_back(p2 - t2);
//                dr.push_back(p3 - t3);
//            }
//            t1 = p1;
//            t2 = p2;
//            t3 = p3;
//        }
//        delta.push_back(dr);
//    }

}

void Compression::pca_encode() {
    Frames frames = clst.get_frames();
    vector<cluster > clusters = clst.get_clusters();
    vector<Matrix > U;
    vector<Matrix > data;
    for (uint32_t i = 0 ; i != clst.get_clusters_num(); ++i) {
        Matrix m(clusters[i].size() * 3, frames.size());
        for (uint32_t j = 0 ; j != clusters[i].size(); ++j) {
            uint32_t p = clusters[i][j];
            for (uint32_t k = 0 ; k != frames.size(); ++k) {
                m(3 * j, k) = frames[k].vertices[p].x;
                m(3 * j + 1, k) = frames[k].vertices[p].y;
                m(3 * j + 2, k) = frames[k].vertices[p].y;
            }
        }
        data.push_back(m);
    }

    for (uint32_t i = 0 ; i != data.size() ; ++i) {
        pca_per_cluster(data[i], 1);
    }
}

void Compression::svd_encode() {
    Frames frames = clst.get_frames();
    vector<cluster > clusters = clst.get_clusters();
    vector<Matrix > U;
    vector<Matrix > data;
    for (uint32_t i = 0 ; i != clst.get_clusters_num(); ++i) {
        Matrix m(clusters[i].size() * 3, frames.size());
        for (uint32_t j = 0 ; j != clusters[i].size(); ++j) {
            uint32_t p = clusters[i][j];
            for (uint32_t k = 0 ; k != frames.size(); ++k) {
                m(3 * j, k) = frames[k].vertices[p].x;
                m(3 * j + 1, k) = frames[k].vertices[p].y;
                m(3 * j + 2, k) = frames[k].vertices[p].y;
            }
        }
        data.push_back(m);
    }

    for (uint32_t i = 0 ; i != data.size() ; ++i) {
        svd_per_cluster(data[i], 1);
    }
}

void Compression::pca_per_cluster(Matrix m, int N) {
    Matrix Fea = m, Fea_reduce;
    vector<vector<double>> mdata(1,Fea.mean(1));
    Matrix mu(1,Fea.Getcols());
    mu.setdata(mdata);
    if(Fea.Getrows()>Fea.Getcols()){
        double size = Fea.Getrows();
        Matrix St = Fea.trans()*Fea/size-mu.trans()*mu;
        Matrix V;
        vector<double> ei;
        if(!St.eig(ei,&V)){
            exit(-1);
        }
        V = V.trunc(N,2);
        Fea_reduce = Fea * V;
        U.push_back(V);
        C.push_back(Fea_reduce);
    }
    else{								//For High-Dimensional conditions
        Matrix *t = new Matrix;
        (*t) = ones(Fea.Getrows(),1)*mu;
        Matrix St = (Fea-(*t))*(Fea.trans()-(*t).trans());	//U*U' -> U'*U
        delete t;
        Matrix Q;
        vector<double> ei;
        if(!St.eig(ei,&Q)){
            exit(-1);
        }
        ei.pop_back();
        Matrix I = diag(ei);
        I = I.m_sqrt();
        I = I.invert();
        Q = Q.trunc(Q.Getcols()-1,2);
        Matrix V = Fea.trans()*Q*I;
        V = V.trunc(N,2);
        Fea_reduce = Fea * V;
        U.push_back(V);
        C.push_back(Fea_reduce);
    }
}

void Compression::svd_per_cluster(Matrix m, int N) {
    Matrix t = m * (m.trans());
    Matrix U, m_reduce;
    vector<double> ei;
    t.eig(ei, &U);
    U = U.trunc(N, 1);
    m_reduce = U * m;
    this->U.push_back(U);
    this->C.push_back(m_reduce);
}


void Compression::quantitzation(int pca_bits, int seed_bits, int delta_bits) {
    double Umax = -MAXFLOAT, Umin = MAXFLOAT, Cmax = -MAXFLOAT, Cmin = MAXFLOAT;
    bithead.pca_bits = pca_bits;
    bithead.delta_bits = delta_bits;
    bithead.seed_bits = seed_bits;
    for (uint32_t i = 0 ; i != U.size(); ++i) {
        double t, p;
        t = U[i].max_element();
        p = U[i].min_element();
        if (t > Umax) Umax = t;
        if (p < Umin) Umin = p;

        U[i].quantization(pca_bits, t - p);
        UL.push_back(t - p);
        Urow_array.push_back(U[i].Getrows());
        Ucol_array.push_back(U[i].Getcols());

        t = C[i].max_element();
        p = C[i].min_element();

        if (t > Cmax) Cmax = t;
        if (t < Cmin) Cmin = p;

        C[i].quantization(pca_bits, t - p);
        CL.push_back(t - p);
        Crow_array.push_back(C[i].Getrows());
        Ccol_array.push_back(C[i].Getcols());
    }

    double d_max = -MAXFLOAT, d_min = MAXFLOAT;
    for (uint32_t i = 0 ; i != seed_triangle_points.size() ; ++i) {
        if (seed_triangle_points[i] > d_max) d_max = seed_triangle_points[i];
        if (seed_triangle_points[i] < d_min) d_min = seed_triangle_points[i];
    }
    for (uint32_t i = 0 ; i != seed_triangle_points.size() ; ++i)
        seed_triangle_points[i] = (seed_triangle_points[i] / (d_max - d_min)) * pow(2, seed_bits);
    bithead.seed_length = d_max - d_min;

    bithead.delta_length = delta.max_element() - delta.min_element();
    delta.quantization(bithead.delta_bits, bithead.delta_length);
}

void Compression::write_bitstream() {
    bithead.seed_num = clst.get_clusters_num();
    bithead.frame_num = clst.get_frames().size();

    in_stream.resize_capacity(1024);
    //write bitstream head
    in_stream.write_bytes(&bithead, sizeof(bithead));
    cout << "write bithead, " << sizeof(bithead) <<" bytes." << endl;

    //write seeds
    for (uint32_t i = 0 ; i != seed_triangle_points.size(); ++i) {
        int32_t t = (int32_t)round(seed_triangle_points[i]);
        in_stream.write_bits(&t, bithead.seed_bits);
    }

    //write delta
    for (uint32_t i = 0 ; i != delta.Getrows(); ++i) {
        for (uint32_t j = 0 ; j != delta.Getcols(); ++j) {
            int32_t t = (int32_t)round(delta(i,j));
            in_stream.write_bits(&t, bithead.delta_bits);
        }
    }

    //write U/C information
    in_stream.write_bytes(UL.data(), UL.size() * sizeof(double));
    in_stream.write_bytes(CL.data(), UL.size() * sizeof(double));
    in_stream.write_bytes(Urow_array.data(), Urow_array.size() * sizeof(int));
    in_stream.write_bytes(Ucol_array.data(), Ucol_array.size() * sizeof(int));
    in_stream.write_bytes(Crow_array.data(), Crow_array.size() * sizeof(int));
    in_stream.write_bytes(Ccol_array.data(), Ccol_array.size() * sizeof(int));

    //write U
    for (uint32_t i = 0 ; i != U.size(); ++i)
        for (uint32_t m = 0 ; m != U[i].Getrows(); ++m)
            for (uint32_t n = 0; n != U[i].Getcols(); ++n) {
                int32_t t = (int32_t)round(U[i](m,n));
                in_stream.write_bits(&t, bithead.pca_bits);
            }

    // write C
    for (uint32_t i = 0 ; i != C.size(); ++i)
        for (uint32_t m = 0 ; m != C[i].Getrows(); ++m)
            for (uint32_t n = 0; n != C[i].Getcols(); ++n) {
                int32_t t = (int32_t)round(C[i](m,n));
                in_stream.write_bits(&t, bithead.pca_bits);
            }
}

double Compression::mse_error(Matrix p, Matrix q) {
    assert(p.Getcols() == q.Getcols());
    assert(p.Getrows() == q.Getrows());
    double error = 0.0;

    for (uint32_t i = 0 ; i != p.Getrows(); ++i)
        for (uint32_t j = 0 ; j != p.Getcols(); ++j)
            error += (p(i, j) - q(i, j)) * (p(i, j) - q(i, j));

    return error / (p.Getcols() * p.Getrows());
}

Matrix Compression::pca_decode(Matrix C, Matrix V) {
    return C * (V.trans());
}

Matrix Compression::svd_decode(Matrix C, Matrix U) {
    return (U.trans() * std::move(C));
}

void Compression::adpac_encoder() {
    entropy_coder.encode(&in_stream, &out_stream);
}

void Compression::rate_distortion() {

}

void Compression::rd_optimization() {

}




