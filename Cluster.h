//
// Created by haines on 3/26/19.
//

#ifndef ANIMATION_COMPRESSION_CLUSTER_H
#define ANIMATION_COMPRESSION_CLUSTER_H

#include <iostream>
#include <cstring>
#include <vector>
#include <map>
#include <math.h>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
#include "WavefrontObjFile.h"


using namespace std;

using Frame = Obj<float_t >;
using Frames = vector<Frame >;
using Point = Point3D<float_t >;
using Matrix33 = Eigen::Matrix3f;
using Vector3 = Point;
using cluster = vector<uint32_t >;

struct LC {
    vector<uint32_t > vertices;
    uint32_t face_idx;
    Matrix33 affline_matrix, restore_matrix;
    Point o;
    Vector3 ix, iy, iz;

    Point world2local(const Point& p);
    Point local2world(const Point& p);

};

struct LCF {
    uint32_t frame_idx;
    vector<LC > lcs;
};


class Cluster {
    Frames obj_frames;
    int num_clusters;
    vector<uint32_t > seeds;
    vector<TriangleFace > seed_triangles;
    vector<LCF > lcfs;
    vector<cluster > clusters;

    uint32_t num_vertices, num_faces;


public:
    Cluster () {}
    Cluster (Frames obj_frames);
    Cluster (Frames obj_frames, int num_clusers);

    void set_frames(Frames obj_frames);
    void set_clusters_num(int N);

    void segmentation();

    const Frames & get_frames() {return obj_frames;}
    const vector<TriangleFace>& get_seed_triangles() {return seed_triangles; }
    int get_clusters_num() {return num_clusters; }
    vector<cluster > get_clusters() {return clusters;}
    vector<LCF > get_lcfs() {return lcfs;}


private:
    void vertices_cluster();
    uint32_t find_face(uint32_t vertex_idx, uint32_t frame_idx);
    LC calc_lcf(uint32_t vertex_idx, uint32_t frame_idx, bool* is_used);
    LC calc_lcf(uint32_t face_idx, uint32_t frame_idx);
    float_t calc_theta(uint32_t vertex_idx, uint32_t lc_idx);

    void select_seeds(uint32_t frame_idx = 0);

    void segmentation_per_frame(const Frame& frame);

    void lcf_construction();
    void lcf_construction_first_frame(uint32_t frame_idx = 0);
    void lcf_construction_per_frame(uint32_t frame_idx, uint32_t prior_frame_idx);

};


#endif //ANIMATION_COMPRESSION_CLUSTER_H
