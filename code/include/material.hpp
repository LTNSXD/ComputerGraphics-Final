#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include <iostream>

#define DIFF    0  // 漫反射
#define SPEC    1  // 镜面反射
#define REFR    2  // 折射

class Material {
public:

    explicit Material(const Vector3f &a_color, 
                      const Vector3f &d_color, 
                      const Vector3f &s_color = Vector3f::ZERO, 
                      int refl_t = DIFF, double refr_index = 1.0 /* 折射率 */) :
            ambientColor(a_color), diffuseColor(d_color), specularColor(s_color), reflectionType(refl_t), refractiveIndex(refr_index) {

    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }

    virtual Vector3f getAmbientColor() const {
        return ambientColor;
    }

    virtual int getReflectionType() const {
        return reflectionType;
    }

    virtual double getRefractiveIndex() const {
        return refractiveIndex;
    }

protected:
    int reflectionType;
    double refractiveIndex;
    Vector3f ambientColor;
    Vector3f diffuseColor;
    Vector3f specularColor;
};


#endif // MATERIAL_H
