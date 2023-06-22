#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include "utils.hpp"
#include <tuple>
#include <cmath>

#define BEGIN       0.0
#define END         1.0
#define MAX_ITER    64
#define EPS         1e-4

using std::min;
using std::max;

class RevSurface : public Object3D {

    Curve *pCurve;

public:
    RevSurface(Curve *pCurve, Material* material) : pCurve(pCurve), Object3D(material) {
        // Check flat.
        Vector3f ld = -vecInf, ru = vecInf;
        for (const auto &cp : pCurve->getControls()) {
            if (cp.z() != 0.0) {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
            Vector3f v(abs(cp.x()), cp.y(), abs(cp.z()));
            ld[0] = min(ld[0], -v[0]), ld[1] = min(ld[1], v[1]), ld[2] = min(ld[2], -v[0]);
            ru[0] = max(ru[0], v[0]), ru[1] = max(ru[1], v[1]), ru[2] = max(ru[2], v[0]);
        }
        this->setAABB(ld, ru);
    }

    ~RevSurface() override {
        delete pCurve;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        // implement this for the ray-tracing routine using G-N iteration.
        if (!this->aabb.intersect(r)) return false;
        if (r.getDirection().y() != 0) { return Intersect(r, h, tmin); }
        else { return SpecialIntersect(r, h, tmin); }
    }

    // bool SpecialIntersect(const Ray &r, Hit &h, double tmin) {
    //     bool res = false;
    //     const Vector3f &o = r.getOrigin(), &d = r.getDirection(), &rd = r.getRdirection();
    //     double d_x_y = d.x() * rd.y(), d_z_y = d.z() * rd.y();
    //     for (int i = 0; i < 128; ++i) {
    //         double t = BEGIN + (END - BEGIN) * (double)i / 128;
    //         for (int j = 0; j < MAX_ITER; ++j) {
    //             const auto &comb = pCurve->pointAtParam(t);
    //             const Vector3f &point = comb.V;
    //             const Vector3f &tangent = comb.T;
    //             double fd = point.y() - o.y();
    //             if (fabs(fd) < EPS) {
    //                 double a = 1.0; // d.x() * d.x() + d.z() * d.z(), for d.y = 0
    //                 double half_b = d.x() * o.x() + d.z() * o.z();
    //                 double c = o.x() * o.x() + o.z() * o.z() - point.x() * point.x();
    //                 double r_t = (-half_b - sqrt(half_b * half_b - c)); // a = 1
    //                 if (tmin < r_t && r_t < h.getT()) {
    //                     double cos_theta = point.x() == 0 ? 0 : (r_t * d.x() + o.x()) / point.x();
    //                     double sin_theta = point.x() == 0 ? 1 : (r_t * d.z() + o.z()) / point.x();
    //                     Vector3f n = Vector3f::cross(t, -Vector3f::FORWARD);
    //                     Vector3f n_rot = Vector3f(
    //                         cos_theta * n.x() - sin_theta * n.z(),
    //                         n.y(),
    //                         sin_theta * n.x() + cos_theta * n.z()
    //                     );
    //                     bool into = true;
    //                     if (Vector3f::dot(n_rot, d) > 0) {
    //                         n_rot = -n_rot;
    //                         into = false;
    //                     }
    //                     h.set(r_t, material, n_rot, into);
    //                     res = true;
    //                 }
    //                 break;
    //             }
    //             double derive = tangent.y();
    //             t -= fd / (32 * derive);
    //             if (t < 0 || t > 1) break;
    //         }
    //     }
    //     return res;
    // }

    // bool Intersect(const Ray &r, Hit &h, float tmin) {
    //     // dy != 0
    //     // Equation: F(t) = (x_cos) ** 2 + (x_sin) ** 2 - x(t) ** 2 = 0
    //     // Tangent: F'(t) = 2 * (x_cos) * (dx / dy) * y'(t) + 
    //     //                  2 * (x_sin) * (dz / dy) * y'(t) -
    //     //                  2 * x(t) * x'(t)
    //     bool res = false;
    //     const Vector3f& o = r.getOrigin(), &d = r.getDirection(), &rd = r.getRdirection();
    //     for (int i = 0; i < 128; ++i) {
    //         double t = BEGIN + (END - BEGIN) * (double)i / 128;
    //         for (int j = 0; j < MAX_ITER; ++j) {
    //             const auto &comb = pCurve->pointAtParam(t);
    //             const Vector3f &point = comb.V;
    //             const Vector3f &tangent = comb.T;
    //             double x_cos = (point.y() - o.y()) * d.x() * rd.y() + o.x(), x_sin = (point.y() - o.y()) * d.z() * rd.y() + o.z(), x = point.x();
    //             double fd = x_cos * x_cos + x_sin * x_sin - x * x;
    //             if (fabs(fd) < EPS) {
    //                 double tr = std::max((x_cos - o.x()) * rd.x(), (x_sin - o.z()) * rd.z());
    //                 if (tmin < tr && tr < h.getT()) {
    //                     double cos_theta = (d.x() * tr + o.x()) / x, sin_theta = (d.z() * tr + o.z()) / x;
    //                     Vector3f Pt(tangent.x() * cos_theta, tangent.y(), tangent.x() * sin_theta), Ptheta(-x * sin_theta, 0, x * cos_theta);
    //                     Vector3f normal = Vector3f::cross(Pt, Ptheta).normalized();
    //                     if (Vector3f::dot(normal, d) > 0) { normal = -normal; }
    //                     h.set(tr, material, normal);
    //                     res = true;
    //                     break;
    //                 }
    //             }
    //             double dif = 2 * (x_cos) * (d.x() / d.y()) * tangent.y() + 
    //                         2 * (x_sin) * (d.z() / d.y()) * tangent.y() - 
    //                         2 * x * tangent.x();
    //             t -= fd / dif;
    //             if (t < 0 || t > 1) { break; }
    //         }
    //     }
    //     return res;
    // }

    bool SpecialIntersect(const Ray &r, Hit &h, double tmin) {
        // Intersect when d.y() == 0
        // equation: y - o.y = 0
        // derive: y'
        // (r_t * d_x + o_x) ** 2 + (r_t * d_z + o_z) ** 2 = x_t ** 2
        bool ret = false;
        const Vector3f &o = r.getOrigin(), &d = r.getDirection(), &rd = r.getRdirection();
        double start = pCurve->start(), end = pCurve->end();
        for (int i = 0; i < 128; ++i) {
            double cur_place = start + (end - start) * double(i) / 128;
            for (int j = 0; j < MAX_ITER; ++j) {
                auto point = pCurve->pointAtParam(cur_place);
                const Vector3f &v = point.V, &t = point.T;
                double solution = v.y() - o.y();
                if (solution > -EPS && solution < EPS) {
                    // solve the r_t equation
                    double a = 1.0; // d.x() * d.x() + d.z() * d.z(), for d.y = 0
                    double half_b = d.x() * o.x() + d.z() * o.z();
                    double c = o.x() * o.x() + o.z() * o.z() - v.x() * v.x();
                    double r_t = (-half_b - sqrt(half_b * half_b - c)); // a = 1
                    if (tmin < r_t && r_t < h.getT() - EPS) {
                        double cos_theta = v.x() == 0 ? 0 : (r_t * d.x() + o.x()) / v.x();
                        double sin_theta = v.x() == 0 ? 1 : (r_t * d.z() + o.z()) / v.x();
                        Vector3f n = Vector3f::cross(t, -Vector3f::FORWARD);
                        Vector3f n_rot = Vector3f(
                            cos_theta * n.x() - sin_theta * n.z(),
                            n.y(),
                            sin_theta * n.x() + cos_theta * n.z()
                        );
                        // assume that normal at x < 0 all to x < 0
                        bool into = true;
                        if (Vector3f::dot(n_rot, d) > 0) {
                            n_rot = -n_rot;
                            into = false;
                        }
                        h.set(r_t, material, n_rot, into);
                        ret = true;
                    }
                    break;
                } // solved
                double derive = t.y();
                cur_place -= solution / (32 * derive);
                if (isnan(cur_place) || cur_place < 0 || cur_place > 1) break;
            }
        }
        return ret;
    }

    bool Intersect(const Ray &r, Hit &h, double tmin) {
        // Intersect when d.y() != 0
        // equation:
        // ((y - o.y) * d.x / d.y + o.x) ** 2 + ((y - o.y) * d.z / d.y + o.z) ** 2 - x ** 2 = 0
        // derive:
        // 2((y - o.y) * d.x / d.y + o.x) * (d.x / d.y) * y' +
        // 2((y - o.y) * d.z / d.y + o.z) * (d.z / d.y) * y' -
        // 2xx'
        bool ret = false;
        const Vector3f &o = r.getOrigin(), &d = r.getDirection(), &rd = r.getRdirection();
        double d_x_y = d.x() * rd.y(), d_z_y = d.z() * rd.y();
        double start = pCurve->start(), end = pCurve->end();
        for (int i = 0; i < 128; ++i) {
            double cur_place = start + (end - start) * double(i) / 128;
            for (int j = 0; j < MAX_ITER; ++j) {
                auto point = pCurve->pointAtParam(cur_place);
                const Vector3f &v = point.V, &t = point.T;
                double x_cos_theta = (v.y() - o.y()) * d_x_y + o.x();
                double x_sin_theta = (v.y() - o.y()) * d_z_y + o.z();
                double solution = x_cos_theta * x_cos_theta + x_sin_theta * x_sin_theta - v.x() * v.x();
                if (solution > -EPS && solution < EPS) {
                    double r_t = max((x_cos_theta - o.x()) * rd.x(), (x_sin_theta - o.z()) * rd.z()); // guarantee no 0 when d.x() == 0 (on Y axis)
                    if (tmin < r_t && r_t < h.getT() - EPS) {
                        double cos_theta = v.x() == 0 ? 0 : x_cos_theta / v.x();
                        double sin_theta = v.x() == 0 ? 1 : x_sin_theta / v.x();
                        Vector3f n = Vector3f::cross(t, -Vector3f::FORWARD);
                        Vector3f n_rot = Vector3f(
                            cos_theta * n.x() - sin_theta * n.z(),
                            n.y(),
                            sin_theta * n.x() + cos_theta * n.z()
                        );
                        // assume that normal at x < 0 all to x < 0
                        bool into = true;
                        if (Vector3f::dot(n_rot, d) > 0) {
                            n_rot = -n_rot;
                            into = false;
                        }
                        h.set(r_t, material, n_rot, into);
                        ret = true;
                    }
                    break;
                } // solved
                double derive = 2 * x_cos_theta * d_x_y * t.y() + 
                                2 * x_sin_theta * d_z_y * t.y() - 2 * v.x() * t.x();
                cur_place -= solution / (32 * derive);
                if (isnan(cur_place) || cur_place < 0 || cur_place > 1) break;
            }
        }
        return ret;
    }
};

#endif //REVSURFACE_HPP
