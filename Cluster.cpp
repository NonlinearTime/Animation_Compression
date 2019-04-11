#include <utility>


//
// Created by haines on 3/26/19.
//

#include "Cluster.h"
#include <utility>

using namespace std;

Cluster::Cluster(Frames obj_frames): obj_frames(std::move(obj_frames)){
    num_vertices = obj_frames[0].num_vertices;
}

void Cluster::set_frames(Frames obj_frames) {
    this->obj_frames = std::move(obj_frames);
    num_vertices = this->obj_frames[0].num_vertices;
}

void Cluster::set_clusters_num(int num_clusters) {
    this->num_clusters = num_clusters;
    clusters.resize(num_clusters);
}

Cluster::Cluster(Frames obj_frames, int num_clusters): obj_frames(move(obj_frames)), num_clusters(num_clusters) {
    clusters.resize(num_clusters);
    num_vertices = obj_frames[0].num_vertices;
}

void Cluster::segmentation_per_frame(const Frame& frame) {

}

uint32_t Cluster::find_face(uint32_t vertex_idx, uint32_t frame_idx) {
    Frame& frame = obj_frames[frame_idx];
    for (uint32_t i = 0 ; i != frame.faces.size(); ++i) {
        if (frame.faces[i].v[0] == vertex_idx) return i;
    }
    return 0;
}

LC Cluster::calc_lcf(uint32_t vertex_idx,  uint32_t frame_idx, bool* is_used) {
    LC lc;
    uint32_t p1, p2, p3;
    Frame& frame = obj_frames[frame_idx];
    for (uint32_t i = 0; i != frame.faces.size(); ++i) {
        if (is_used[i]) continue;
        p1 = frame.faces[i].v[0];
        p2 = frame.faces[i].v[1];
        p3 = frame.faces[i].v[2];
        if (p1 == vertex_idx || p2 == vertex_idx || p3 == vertex_idx) {
            is_used[i] = true;
            lc.face_idx = i;
            lc.o = (frame.vertices[p1] + frame.vertices[p2]) / 2;

            Vector3 p12 = frame.vertices[p2] - frame.vertices[p1];
            Vector3 p13 = frame.vertices[p3] - frame.vertices[p1];
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

            break;
        }
    }

    return lc;
}

LC Cluster::calc_lcf(uint32_t face_idx, uint32_t frame_idx) {
    LC lc;
    uint32_t p1, p2, p3;
    Frame& frame = obj_frames[frame_idx];
    p1 = frame.faces[face_idx].v[0];
    p2 = frame.faces[face_idx].v[1];
    p3 = frame.faces[face_idx].v[2];

    lc.face_idx = face_idx;
    lc.o = (frame.vertices[p1] + frame.vertices[p2]) / 2;

    Vector3 p12 = frame.vertices[p2] - frame.vertices[p1];
    Vector3 p13 = frame.vertices[p3] - frame.vertices[p1];
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

    return lc;
}

void Cluster::vertices_cluster() {
    float_t min_theta = MAXFLOAT;
    float_t theta;
    uint32_t cluster_idx = 0;
    bool is_visited[num_vertices];
    memset(is_visited, 0, sizeof(bool) * num_vertices);
    uint32_t sum = 0;

    for (uint32_t j = 0 ; j != num_clusters; ++j) {
        TriangleFace& f = obj_frames[0].faces[lcfs[0].lcs[j].face_idx];
        for (uint32_t i = 0; i != 3; ++i) {
            if (!is_visited[f.v[i]]) {
                sum += 1;
                is_visited[f.v[i]] = true;
            }
        }
    }

    cout << "seed points num: " << sum << endl;

    for (uint32_t i = 0 ; i != num_vertices; ++i) {
        min_theta = MAXFLOAT;
        for (uint32_t j = 0 ; j != num_clusters; ++j) {
            theta = calc_theta(i, j);
            if (theta < min_theta) {
                cluster_idx = j;
                min_theta = theta;
            }
        }
//        TriangleFace& face = obj_frames[0].faces[lcfs[0].lcs[cluster_idx].face_idx];
//        if (i != face.v[0] && i != face.v[1] && i != face.v[2])
//            clusters[cluster_idx].push_back(i);
        if (!is_visited[i]) clusters[cluster_idx].push_back(i);
    }

    sum = 0;
    for (uint32_t i = 0; i != num_clusters; ++i) {
        seed_triangles.push_back(obj_frames[0].faces[lcfs[0].lcs[i].face_idx]);
        sum += clusters[i].size();
        cout << clusters[i].size() << endl;
    }
    cout << "Total clusters vetices num: " << sum << endl;
}

float_t Cluster::calc_theta(uint32_t vertex_idx, uint32_t lc_idx) {
    vector<float_t > world_coord;
    float_t theta = 0;
    Point prior_coord, current_local_coord, current_world_coord;
    prior_coord = lcfs[0].lcs[lc_idx].world2local(obj_frames[0].vertices[vertex_idx]);

    for (uint32_t i = 1 ; i != obj_frames.size(); ++i) {
        current_world_coord = obj_frames[i].vertices[vertex_idx];
        current_local_coord = lcfs[i].lcs[lc_idx].world2local(current_world_coord);
        theta += vec_length2(current_local_coord - prior_coord);
        prior_coord = current_local_coord;
    }

    return theta;
}

void Cluster::reduce_vertices(uint32_t vertex_idx, uint32_t frame_idx, bool *is_used) {
    Frame& frame = obj_frames[frame_idx];
//    for (uint32_t i = 0 ; i != frame.faces.size(); ++i) {
//        if (frame.faces[i].v[0] == vertex_idx) {
//            is_used[frame.faces[i].v[1]] = true;
//            is_used[frame.faces[i].v[2]] = true;
//        }
//        if (frame.faces[i].v[1] == vertex_idx) {
//            is_used[frame.faces[i].v[0]] = true;
//            is_used[frame.faces[i].v[2]] = true;
//        }
//        if (frame.faces[i].v[2] == vertex_idx) {
//            is_used[frame.faces[i].v[0]] = true;
//            is_used[frame.faces[i].v[1]] = true;
//        }
//    }

    for (uint32_t i = 0; i != frame.topo[vertex_idx].size(); ++i) {
        uint32_t v = frame.topo[vertex_idx][i];
        is_used[v] = true;
        for (uint32_t j = 0; j != frame.topo[v].size(); ++j) {
            is_used[frame.topo[v][j]] = true;
        }
    }
}


void Cluster::select_seeds(uint32_t frame_idx) {
    Frame& frame = obj_frames[frame_idx];
    Point3D<float_t > frame_center = frame.bounding_sphere_c;
    uint64_t num_vertices = frame.vertices.size();
    bool is_visited[num_vertices];
    memset(is_visited, 0, sizeof(bool) * num_vertices);

    //find seeds
    float_t max_dist = -1;
    uint32_t idx = 0;

    for (uint32_t i = 0 ; i != num_vertices; ++i) {
        if (dis(frame.vertices[i], frame_center) > max_dist) {
            idx = i;
            max_dist = dis(frame.vertices[i], frame_center);
        }
    }

    is_visited[idx] = true;
    reduce_vertices(idx, frame_idx, is_visited);
    seeds.push_back(idx);

    float_t min_dist = MAXFLOAT;
    max_dist = -1;

    int n = num_clusters - 1;
    while (n--) {
        max_dist = -1;
        for (uint32_t i = 0; i != num_vertices; ++i) {
            if (is_visited[i]) continue;
            min_dist = MAXFLOAT;
            for (uint32_t j = 0; j != seeds.size(); ++j) {
                float_t seed_dist = dis(frame.vertices[seeds[j]], frame.vertices[i]);
                if (seed_dist < min_dist) {
                    min_dist = seed_dist;
                }
            }
            if (min_dist > max_dist) {
                idx = i;
                max_dist = min_dist;
            }
        }
        is_visited[idx] = true;
        reduce_vertices(idx, frame_idx, is_visited);
        seeds.push_back(idx);
    }

    cout << "seeds number: " << seeds.size() << endl;
}

void Cluster::lcf_construction_first_frame(uint32_t frame_idx) {
    //find affline matrix, size: 3 * 3
    LCF lcf;
    lcf.frame_idx = frame_idx;
    bool is_used[obj_frames[frame_idx].faces.size()];
    memset(is_used, 0, sizeof(bool) * obj_frames[frame_idx].faces.size());
    for (uint32_t i = 0 ; i != seeds.size(); ++i) {
        lcf.lcs.push_back(calc_lcf(seeds[i], frame_idx, is_used));
    }
    lcfs.push_back(lcf);
}

void Cluster::lcf_construction_per_frame(uint32_t frame_idx, uint32_t prior_frame_idx) {
    assert(!lcfs[prior_frame_idx].lcs.empty());
    uint32_t face_idx ;
    LCF lcf;
    lcf.frame_idx = frame_idx;
    for (uint32_t i = 0 ; i != lcfs[prior_frame_idx].lcs.size(); ++i) {
        face_idx = lcfs[prior_frame_idx].lcs[i].face_idx;
        lcf.lcs.push_back(calc_lcf(face_idx, frame_idx));
    }
    lcfs.push_back(lcf);
}

void Cluster::lcf_construction() {
    lcf_construction_first_frame(0);
    for (uint32_t i = 1 ; i != obj_frames.size() ; ++i) {
        lcf_construction_per_frame(i, i - 1);
    }
}

void Cluster::segmentation() {
    select_seeds();
    lcf_construction();
    vertices_cluster();
}


Point LC::world2local(const Point& p) const {
    Point l, pp;
    pp = p;
    pp -= o;
    l.x = affline_matrix(0,0) * pp.x + affline_matrix(0,1) * pp.y + affline_matrix(0,2) * pp.z;
    l.y = affline_matrix(1,0) * pp.x + affline_matrix(1,1) * pp.y + affline_matrix(1,2) * pp.z;
    l.z = affline_matrix(2,0) * pp.x + affline_matrix(2,1) * pp.y + affline_matrix(2,2) * pp.z;
    return l;
}

Point LC::local2world(const Point& p) const {
    Point l;
    l.x = restore_matrix(0,0) * p.x + restore_matrix(0,1) * p.y + restore_matrix(0,2) * p.z;
    l.y = restore_matrix(1,0) * p.x + restore_matrix(1,1) * p.y + restore_matrix(1,2) * p.z;
    l.z = restore_matrix(2,0) * p.x + restore_matrix(2,1) * p.y + restore_matrix(2,2) * p.z;
    l += o;
    return l;
}
