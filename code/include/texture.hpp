#ifndef TEXTURE_H
#define TEXTURE_H

#include <iostream>
#include <vecmath.h>
#include "utils.hpp"
#include "image.hpp"

class Texture {

    Image *image;

public: 
    Texture(const char *filename) {
        image = Image::LoadBMP(filename);
    }

    Vector3f getPixel(const Vector2f &pos) {
        int w = image->Width(), h = image->Height();
        double x = pos.x(), y = pos.y();
        x -= (int)x; y -= (int)y;
        x = x > 0. ? x : .99 + x;
        y = y > 0. ? y : .99 + y;
        return image->GetPixel(x * w, y * h);
    }
};



#endif