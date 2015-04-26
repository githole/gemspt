#ifndef	_VEC_H_
#define	_VEC_H_

#include <cmath>

#include "constant.h"

namespace gemspt {
    
// ベクトル演算用クラス
struct Vec {
    double x, y, z;

    Vec(const double x = 0, const double y = 0, const double z = 0) : x(x), y(y), z(z) {}
    
    inline Vec operator+(const Vec &v) const {
        return Vec(x + v.x, y + v.y, z + v.z);
    }
    inline Vec operator-(const Vec &v) const {
        return Vec(x - v.x, y - v.y, z - v.z);
    }
    inline Vec operator*(const double a) const {
        return Vec(x * a, y * a, z * a);
    }
    inline Vec operator/(const double a) const {
        return Vec(x / a, y / a, z / a);
    }
    inline const double length_squared() const { 
        return x*x + y*y + z*z; 
    }
    inline Vec operator-() const {
        return Vec(-x, -y, -z);
    }
    inline const double length() const { 
        return sqrt(length_squared()); 
    }
};

inline Vec operator*(const double a, const Vec &v) { 
    return v * a; 
}

inline Vec normalize(const Vec &v) {
    return v * (1.0 / v.length()); 
}

inline Vec multiply(const Vec &v1, const Vec &v2) {
    return Vec(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

inline double dot(const Vec &v1, const Vec &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline Vec cross(const Vec &v1, const Vec &v2) {
    return Vec(
        (v1.y * v2.z) - (v1.z * v2.y),
        (v1.z * v2.x) - (v1.x * v2.z),
        (v1.x * v2.y) - (v1.y * v2.x));
}

inline Vec reflect(const Vec &in, const Vec &normal) {
    return normalize(in - normal * 2.0 * dot(normal, in));
}

// 正規直交基底を作る
inline void createOrthoNormalBasis(const Vec &normal, Vec *tangent, Vec *binormal) {
    if (abs(normal.x) > abs(normal.y))
    *tangent = normalize(cross(Vec(0, 1, 0), normal));
    else
    *tangent = normalize(cross(Vec(1, 0, 0), normal));
    *binormal = normalize(cross(normal, *tangent));
}

};

#endif
