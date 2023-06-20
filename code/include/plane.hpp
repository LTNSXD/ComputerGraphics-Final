#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane() = delete;

    Plane(const Vector3f &normal, float d, Material *m) : Object3D(m), normal(normal.normalized()), d(d) {
        
    }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        Vector3f R0 = r.getOrigin(), Rd = r.getDirection();
        float t = - ( - d + Vector3f::dot(normal, R0)) / (Vector3f::dot(normal, Rd));
        if (t >= 0 && t > tmin && t < h.getT()) {
            h.set(t, material, normal);
            return true;
        }
        else {
            return false;
        }
    }

protected:
    Vector3f normal;
    float d;

};

#endif //PLANE_H
		

