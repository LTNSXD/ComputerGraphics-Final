// #ifndef CAMERA_H
// #define CAMERA_H

// #include "ray.hpp"
// #include "utils.hpp"
// #include <vecmath.h>
// #include <float.h>
// #include <cmath>
// #include <iostream>


// class Camera {
// public:
//     Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
//         this->center = center;
//         this->direction = direction.normalized();
//         this->horizontal = Vector3f::cross(this->direction, up).normalized();
//         this->up = Vector3f::cross(this->horizontal, this->direction);
//         this->width = imgW;
//         this->height = imgH;
//     }

//     // Generate rays for each screen-space coordinate
//     virtual Ray generateRay(const Vector2f &point) = 0;
//     virtual ~Camera() = default;

//     int getWidth() const { return width; }
//     int getHeight() const { return height; }

// protected:
//     // Extrinsic parameters
//     Vector3f center;
//     Vector3f direction;
//     Vector3f up;
//     Vector3f horizontal;
//     // Intrinsic parameters
//     int width;
//     int height;
// };

// // TODO: Implement Perspective camera
// // You can add new functions or variables whenever needed.
// class PerspectiveCamera : public Camera {

// public:
//     PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
//             const Vector3f &up, int imgW, int imgH, float angle, float focal = 1.0, float aperture = 0.0) : Camera(center, direction, up, imgW, imgH) {
//         // angle is in radian.
//         this->angle = angle;
//         this->focal = focal;
//         this->aperture = aperture;
//     }

//     Ray generateRay(const Vector2f &point) override {
//         // 
//         float f = height / (2 * (tan(angle / 2)));
//         float cx = width / 2;
//         float cy = height / 2;

//         Vector3f Rc = Vector3f((point[0] - cx) / f, (cy - point[1]) / f, 1).normalized();

//         Vector3f delta = this->aperture * Vector3f((rand01() - .5) * 2, (rand01() - .5) * 2, 0);

//         Matrix3f R = Matrix3f(horizontal, - up, direction);
//         Vector3f Rw = (R * (focal * Rc + delta)).normalized();

//         return Ray(center, Rw);
//     }

// private:
//     float angle;
//     float focal, aperture;
// };

// #endif //CAMERA_H

#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include "utils.hpp"
#include <vecmath.h>
#include <float.h>
#include <cmath>

class Camera {
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up).normalized();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// TODO: Implement Perspective camera
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {
    float fxy;
public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle) : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        this->fxy = imgH / (2 * tan(angle / 2));
    }

    Ray generateRay(const Vector2f &point) override {
        // 
        Vector3f drc(point.x() - width/2, height/2 - point.y(), fxy);
        drc.normalize();
        Vector3f drw = drc.x() * horizontal - drc.y() * up + drc.z() * direction;
        drw.normalize();
        return Ray(this->center, drw);
    }
};

class LensCamera : public Camera {
    float fxy;
    float distance;//成像平面到其的距离
    float radius;//相机光圈半径
public:
    LensCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle, float distance, float radius) : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        this->fxy = imgH / (2 * tan(angle / 2) * (distance + 1));
        this->distance = distance;//在真实空间的距离
        this->radius = radius;//在真实空间的半径
    }

    Ray generateRay(const Vector2f &point) override {
        // 

        Vector3f drc(point.x() - width/2, height/2 - point.y(), fxy * (distance + 1));
        drc.normalize();
        Vector3f drw = drc.x() * horizontal - drc.y() * up + drc.z() * direction;
        drw.normalize();
        // Vector3f drc(point.x() - width/2, height/2 - point.y(), fxy);
        // drc.normalize();
        float theta = 2 * M_PI * rand01();
        float r = rand01() * radius;
        Vector3f point2 = r * sin(theta) * up +  r * cos(theta) * horizontal;
        Vector3f r2 = drw * ((distance + 1) / (Vector3f::dot(drw, direction))) - point2;
        r2.normalize();
        return Ray(this->center + point2, r2);
        // Vector3f r2 = drc * (distance / fxy) - point2;
        // r2.normalize();
        // Vector3f drw = r2.x() * horizontal - r2.y() * up + r2.z() * direction;
        // Vector3f tran_point2 = point2.x() * horizontal - point2.y() * up;
        // drw.normalize();
        // return Ray(this->center + tran_point2, drw);
    }
};

#endif //CAMERA_H
