#ifndef	_SPHERE_H_
#define	_SPHERE_H_

#include <cmath>

#include "vec.h"
#include "ray.h"
#include "constant.h"
#include "hitpoint.h"

namespace gemspt {

// 球の幾何学的な情報を持つ
class Sphere {
private:
    double radius_;
    Vec position_;
public:
    Sphere(const double radius, const Vec &position) :
      radius_(radius), position_(position) {}

    // 入力のrayに対する交差点までの距離を得る。
    // 交差したらtrue,さもなくばfalseを返す。
    inline bool intersect(const Ray &ray, Hitpoint *hitpoint) const {
        // 自己交差の判定用定数。
        const double kEPS = 1e-6; 

        const Vec o_to_p = position_ - ray.org;
        const double b = dot(o_to_p, ray.dir);
        const double c = b * b - dot(o_to_p, o_to_p) + radius_ * radius_;

        if (c < 0.0)
            return false;
        
        const double sqrt_c = sqrt(c);
        const double t1 = b - sqrt_c, t2 = b + sqrt_c;
    
        // 微小距離内だったら交差しないとする（自己交差を避けるため）。
        if (t1 < kEPS && t2 < kEPS)
            return false;

        // 交差するときは二点以上で交差する。（接する場合は一点）
        // 近い方を交差点とする。また、負値の場合はレイの逆方向で交差してるため交差していないとする。
        if (t1 > kEPS) {
            hitpoint->distance = t1;
        } else {
            hitpoint->distance = t2;
        }

        hitpoint->position = ray.org + hitpoint->distance * ray.dir;
        hitpoint->normal   = normalize(hitpoint->position - position_);
        return true;
    }
};

};

#endif
