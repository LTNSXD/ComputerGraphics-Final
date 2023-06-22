#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>

// TODO: Implement functions and add more fields as necessary

class Sphere : public Object3D {
public:
    Sphere() {
        // unit ball at the center
        center = Vector3f(0, 0, 0);
        radius = 1;
    }

    Sphere(const Vector3f &center, float radius, Material *material) : Object3D(material, center - radius * 1.415, center + radius * 1.415), center(center), radius(radius) {

    }

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        Vector3f l = center - r.getOrigin();
        float tp = Vector3f::dot(l, r.getDirection().normalized());

        if ((l.squaredLength() > radius * radius) && (tp < 0)) {
            // tp < 0
            return false;
        }

        float d_square = l.squaredLength() - tp * tp;
        if (d_square > radius * radius) {
            // d > r
            return false;
        }
        float tr = sqrt(radius * radius - d_square);
        float t = (l.squaredLength() > radius * radius) ? (tp - tr) : (tp + tr);
        // Vector3f intersection = r.pointAtParameter(t);
        Vector3f intersection = r.getOrigin() + t * r.getDirection().normalized();
        Vector3f N = ((l.squaredLength() > radius * radius ? 1 : -1) * (intersection - center)).normalized();
        bool into = l.squaredLength() > radius * radius;
        t = t / r.getDirection().length();
        if (t < h.getT() && t > tmin) {
            h.set(t, material, N, into);
            return true;
        }
        else {
            return false;
        }
    }

    Vector3f getCenter() const { return center; }

    float getRadius() const { return radius; }

protected:
    Vector3f center;
    float radius;
};


#endif
