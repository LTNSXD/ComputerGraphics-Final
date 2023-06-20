#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>
using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle: public Object3D {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) : Object3D(m) {
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;
		// 
		normal = Vector3f::cross((b - a), (c - a)).normalized();
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
		// Implement Crammer's Rule to calculate intersection
		// TODO: Make it faster
        Vector3f E1 = vertices[0] - vertices[1];
		Vector3f E2 = vertices[0] - vertices[2];
		Vector3f S = vertices[0] - ray.getOrigin();
		Vector3f Rd = ray.getDirection();

		Matrix3f M0 = Matrix3f(Rd, E1, E2);
		Matrix3f M1 = Matrix3f(S, E1, E2);
		Matrix3f M2 = Matrix3f(Rd, S, E2);
		Matrix3f M3 = Matrix3f(Rd, E1, S);

		// Cramer
		if (M0.determinant() == 0) {
			return false;
		}
		float t = M1.determinant() / M0.determinant();
		float beta = M2.determinant() / M0.determinant();
		float gamma = M3.determinant() / M0.determinant();

		if (0 < t && 0 < beta && 0 < gamma && beta < 1 && gamma < 1 && beta + gamma < 1 && 
			t < hit.getT() && tmin < t) {
			Vector3f N = (Vector3f::dot(normal, ray.getDirection()) < 0) ? normal : - normal;
			hit.set(t, material, N);
			return true;
		}
		else {
			return false;
		}
	}
	Vector3f normal;
	Vector3f vertices[3];
protected:

};

#endif //TRIANGLE_H
