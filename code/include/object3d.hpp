#ifndef OBJECT3D_H
#define OBJECT3D_H

#include "ray.hpp"
#include "hit.hpp"
#include "material.hpp"
#include "aabb.hpp"

// Base class for all 3d entities.
class Object3D {
public:
    Object3D() : material(nullptr) {}

    virtual ~Object3D() = default;

    Object3D(Material *material) {
        this->material = material;
        this->aabb = AABB();
    }

    Object3D(Material *material, Vector3f ld, Vector3f ru) {
        this->material = material;
        this->aabb = AABB(ld, ru);
    }

    AABB getAABB() {
        return aabb;
    }

    void setAABB(Vector3f ld, Vector3f ru) {
        aabb.set(ld, ru);
    }

    Material *getMaterial() { return material; }

    // Intersect Ray with this object. If hit, store information in hit structure.
    virtual bool intersect(const Ray &r, Hit &h, float tmin) = 0;
protected:
    AABB aabb;
    Material *material;
};

#endif

