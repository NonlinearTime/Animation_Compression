//
// Created by haines on 3/25/19.
//

#ifndef ANIMATION_COMPRESSION_WAVEFRONTOBJFILE_H
#define ANIMATION_COMPRESSION_WAVEFRONTOBJFILE_H

#include <iostream>
#include <cstring>
#include <math.h>
#include <vector>
#include <map>
#include <fstream>

using namespace std;


#define inner_product(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z*(v).z)
#define vec_length2(v)   inner_product(v,v)        // vec_length2 = squared length of vector
#define vec_length(v)    sqrt(vec_length2(v))  // vec_length = length of vector
#define dis(u,v)     vec_length(u-v)       // distance = vec_length of difference

template <typename T>
class Point3D {

public:
    T x, y, z;
    Point3D(T x, T y, T z): x(x), y(y), z(z) {}
    Point3D() {x = y = z = 0;}
    ~Point3D() {}

    Point3D operator+ (Point3D m) {
        Point3D r;
        r.x = x + m.x;
        r.y = y + m.y;
        r.z = z + m.z;
        return r;
    };

    Point3D operator- (Point3D m) {
        Point3D r;
        r.x = x - m.x;
        r.y = y - m.y;
        r.z = z - m.z;
        return r;
    };

    Point3D operator* (Point3D m) {
        Point3D r;
        r.x = x * m.x;
        r.y = y * m.y;
        r.z = z * m.z;
        return r;
    }

    friend Point3D operator* (Point3D p, T m) {
        return Point3D(p.x * m, p.y * m, p.z * m);
    }

    friend Point3D operator* (T m, Point3D p) {
        return Point3D(p.x * m, p.y * m, p.z * m);
    }

    Point3D operator/ (T m) {
        Point3D r;
        r.x = x / m;
        r.y = y / m;
        r.z = z / m;
        return r;
    }


    Point3D& operator= (Point3D m) {
        x = m.x;
        y = m.y;
        z = m.z;
        return *this;
    }

    Point3D& operator+= (const Point3D &m) {
        x += m.x;
        y += m.y;
        z += m.z;
        return *this;
    }

    Point3D& operator-= (const Point3D &m) {
        x -= m.x;
        y -= m.y;
        z -= m.z;
        return *this;
    }

    bool operator== (const Point3D &m) {
        return x == m.x && y == m.y && z == m.z;
    }

    bool operator!= (const Point3D &m) {
        return x != m.x || y != m.y || z != m.z;
    }

    void Print() {
        cout << x << " " << y << " " << z << endl;
    }

    friend ostream &operator<<( ostream &output,
                                const Point3D &m )
    {
        output << m.x << " " << m.y << " " << m.z;
        return output;
    }

    friend istream &operator>>( istream  &input, Point3D &m )
    {
        input >> m.x >> m.y >> m.z;
        return input;
    }

    Point3D& cross_product(const Point3D &x, const Point3D &y) {
        Point3D r;
        r.x = x.y * y.z - y.x * x.y;
        r.y = x.z * y.x - x.x * x.z;
        r.z = x.x * y.y - x.y - y.x;
        return r;
    }

    friend Point3D operator^ (const Point3D& x, const Point3D& y) {
        Point3D r;
        r.x = x.y * y.z - y.x * x.y;
        r.y = x.z * y.x - x.x * x.z;
        r.z = x.x * y.y - x.y - y.x;
        return r;
    }


private:

};

struct TriangleFace{
    uint32_t v[3];

    friend ostream &operator<<( ostream &output,
                                const TriangleFace &t )
    {
        output << t.v[0] << " " << t.v[1] << " " << t.v[2];
        return output;
    }
};


template <typename T>
struct Obj {
    vector<Point3D<T>> vertices;
    vector<TriangleFace> faces;
    Point3D<T> bounding_box[2];
    Point3D<T> bounding_sphere_c;
    T bounding_sphere_r;
    uint32_t num_vertices, num_faces;
    vector<vector<uint32_t >> topo;
};

template <typename T>
class WavefrontObjFile {
    Obj<T > obj;
public:
    WavefrontObjFile() {}
    bool load_obj(const string& file_name);
    void save_obj(const string& file_name);
    Obj<T> get_obj();
    void set_obj(Obj<T> obj);
};

template<typename T>
bool WavefrontObjFile<T>::load_obj(const string& file_name) {
    ifstream fin;
    fin.open(file_name, ios::in);
    if (!fin.good()) {
        cout << "Error: load obj " << file_name << " failed!" << endl;
        return false;
    }

    obj.num_faces = 0;
    obj.num_vertices = 0;
    obj.vertices.clear();
    obj.faces.clear();

    char buf[256], str[256];
    T x, y, z;
    uint32_t id1, id2, id3;
    int line = 0;

    while (!fin.getline(buf, 255).eof()) {
        buf[255] = '\0';

        sscanf(buf, "%s ", str);
        if (strcmp(str, "v") == 0) {
            if (sscanf(buf, "v %f %f %f", &x, &y, &z) == 3) {
                Point3D<float> p(x, y, z);
                obj.vertices.push_back(p);
                obj.num_vertices++;
            } else {
                cout << "Error: vertex format error in file " << file_name << ", line " << line << endl;
                return false;
            }
        } else if (strcmp(str, "f") == 0) {
            if (sscanf(buf, "f %u %u %u", &id1, &id2, &id3) == 3) {
                TriangleFace face;
                face.v[0] = id1;
                face.v[1] = id2;
                face.v[2] = id3;
                obj.faces.push_back(face);
                obj.num_faces++;
            } else {
                cout << "Error: face format error in file " << file_name << ", line " << line << endl;
                return false;
            }
        } else if (strcmp(str, "vt") == 0) {

        } else if (strcmp(str, "vn") == 0) {

        } else if (strcmp(str, "vp") == 0) {

        } else if (strcmp(str, "l") == 0) {

        }
        line++;
    }

    T xmin,ymin,zmin,xmax,ymax,zmax;
    int Pxmin,Pxmax,Pymin,Pymax,Pzmin,Pzmax;

    //calculate the bounding sphere
    xmin = xmax = obj.vertices[0].x;
    ymin = ymax = obj.vertices[0].y;
    zmin = zmax = obj.vertices[0].z;
    Pxmin = Pxmax = Pymin = Pymax = Pzmin = Pzmax = 0;

    //calculate the bounding box
    obj.bounding_box[0] = Point3D<T> (obj.vertices[0].x,obj.vertices[0].y,obj.vertices[0].z);
    obj.bounding_box[1] = Point3D<T> (obj.vertices[0].x,obj.vertices[0].y,obj.vertices[0].z);


    for(unsigned int i = 1; i < obj.vertices.size(); i++)
    {
        //update min value
        obj.bounding_box[0].x = min(obj.vertices[i].x,obj.bounding_box[0].x);
        obj.bounding_box[0].y = min(obj.vertices[i].y,obj.bounding_box[0].y);
        obj.bounding_box[0].z = min(obj.vertices[i].z,obj.bounding_box[0].z);

        //update max value
        obj.bounding_box[1].x = max(obj.vertices[i].x,obj.bounding_box[1].x);
        obj.bounding_box[1].y = max(obj.vertices[i].y,obj.bounding_box[1].y);
        obj.bounding_box[1].z = max(obj.vertices[i].z,obj.bounding_box[1].z);

        //update the  x min and max
        if (obj.vertices[i].x < xmin){
            xmin = obj.vertices[i].x;
            Pxmin = i;
        }
        else if(obj.vertices[i].x > xmax){
            xmax = obj.vertices[i].x;
            Pxmax = i;
        }
        //update the y min and max
        if (obj.vertices[i].y < ymin){
            ymin = obj.vertices[i].y;
            Pymin = i;
        }
        else if(obj.vertices[i].y > ymax){
            ymax = obj.vertices[i].y;
            Pymax = i;
        }
        //update the z min and max
        if(obj.vertices[i].z < zmin){
            zmin = obj.vertices[i].z;
            Pzmin = i;
        }
        else if(obj.vertices[i].z > zmax){
            zmax = obj.vertices[i].z;
            Pzmax = i;
        }
    }

    //calculate the bounding sphere
    Point3D<T>  dVx = obj.vertices[Pxmax] - obj.vertices[Pxmin];
    Point3D<T>  dVy = obj.vertices[Pymax] - obj.vertices[Pymin];
    Point3D<T>  dVz = obj.vertices[Pzmax] - obj.vertices[Pzmin];
    T dx2 = vec_length2(dVx);
    T dy2 = vec_length2(dVy);
    T dz2 = vec_length2(dVz);

    Point3D<T>  center;
    T  radius2;
    T  radius;

    if (dx2 >= dy2 && dx2>=dz2) {					// x direction is largest extent
        center = obj.vertices[Pxmin] + (dVx / 2.0);	// Center = midpoint of extremes
        radius2 = vec_length2(obj.vertices[Pxmax] - center);// radius squared
    }
    else if (dy2 >= dx2  && dy2>=dz2){				// y direction is largest extent
        center = obj.vertices[Pymin] + (dVy / 2.0);	// Center = midpoint of extremes
        radius2 = vec_length2(obj.vertices[Pymax] - center);// radius squared
    }
    else{
        center = obj.vertices[Pzmin] + (dVz / 2.0);	// Center = midpoint of extremes
        radius2 = vec_length2(obj.vertices[Pzmax] - center);// radius squared
    }

    radius = sqrt(radius2);

    // now check that all points V[i] are in the ball
    // and if not, expand the ball just enough to include them
    Point3D<T>  dV;
    float dist2,dist;
    for (unsigned int i=0;i<obj.vertices.size();i++)
    {
        dV = obj.vertices[i] - center;
        dist2 = vec_length2(dV);
        if (dist2 <= radius2) // V[i] is inside the ball already
            continue;

        // V[i] not in ball, so expand ball to include it
        dist = sqrt(dist2);
        radius = (radius + dist) / 2.0;         // enlarge radius just enough
        radius2 = radius * radius;
        center = center + ((dist-radius)/dist) * dV;   // shift Center toward V[i]

    }

    obj.bounding_sphere_c = center;
    obj.bounding_sphere_r = radius;

//    cout<< "----------obj file loaded-------------"<<endl;
    cout<<"number of faces:" << obj.faces.size() << " number of vertices:" << obj.vertices.size() << endl;
//    cout<< "obj bounding box: min:("
//        << obj.bounding_box[0].x << "," << obj.bounding_box[0].y << "," << obj.bounding_box[0].z <<") max:("
//        << obj.bounding_box[1].x << "," << obj.bounding_box[1].y << "," << obj.bounding_box[1].z <<")" << endl
//        <<"obj bounding sphere center:"<<obj.bounding_sphere_c.x<<","<<obj.bounding_sphere_c.y<<","<<obj.bounding_sphere_c.z<<endl
//        <<"obj bounding sphere radius:"<<obj.bounding_sphere_r<<endl;
    obj.topo.resize(obj.num_vertices);
    uint32_t p1, p2, p3;
    for (uint32_t i = 0; i != obj.num_faces; ++i) {
        p1 = obj.faces[i].v[0];
        p2 = obj.faces[i].v[1];
        p3 = obj.faces[i].v[2];

        obj.topo[p1].push_back(p2);
        obj.topo[p1].push_back(p3);

        obj.topo[p2].push_back(p1);
        obj.topo[p2].push_back(p3);

        obj.topo[p3].push_back(p1);
        obj.topo[p3].push_back(p2);
    }

    fin.close();
    return false;
}

template<typename T>
void WavefrontObjFile<T>::save_obj(const string& file_name) {
    ofstream fout;
    fout.open(file_name, ios::trunc | ios::out);

    for (uint32_t i = 0 ; i != obj.vertices.size() ; ++i) {
        fout << "v " << obj.vertices[i] << endl;
    }

    for (uint32_t i = 0 ; i != obj.faces.size() ; ++i) {
        fout << "f " << obj.faces[i] << endl;
    }

    fout.close();
}

template<typename T>
Obj<T> WavefrontObjFile<T>::get_obj() {
    return obj;
}

template<typename T>
void WavefrontObjFile<T>::set_obj(Obj<T> obj) {
    this->obj = obj;
}

#endif //ANIMATION_COMPRESSION_WAVEFRONTOBJFILE_H
