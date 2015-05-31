#ifndef _RADIANCE_H_
#define _RADIANCE_H_

#include <algorithm>

#include "ray.h"
#include "scene.h"
#include "sphere.h"
#include "hitpoint.h"
#include "random.h"

namespace gemspt {

// ray方向からの放射輝度を求める
Color radiance(const Ray &ray, Random &random, const int depth) {
    const Color kBackgroundColor = Color(0.0f, 0.0f, 0.0f);
    const int kDepthLimit = 10;
    // 打ち切りチェック
    if (depth >= kDepthLimit)
        return Color();
    
    // シーンと交差判定
    Hitpoint hitpoint;
    const SceneSphere *now_object = intersect_scene(ray, &hitpoint);
    // 交差チェック
    if (now_object == NULL)
        return kBackgroundColor;

    // マテリアル取得
    const Material *now_material = now_object->get_material();
    const Color emission = now_material->emission();
    if (emission.x > 0.0 || emission.y > 0.0 || emission.z > 0.0) {
        // 光源にヒットしたら放射項だけ返して終わる。
        // （今回、光源は反射率0と仮定しているため）
        return emission;
    }
    
    // 次の方向をサンプリング + その方向のBRDF項の値を得る。
    double pdf = -1;
    Color brdf_value;
    const Vec dir_out = now_material->sample(random, ray.dir, hitpoint.normal, &pdf, &brdf_value);

    // cos項。
    const double cost = dot(hitpoint.normal, dir_out);

    // レンダリング方程式をモンテカルロ積分によって再帰的に解く。
    const Color L = multiply(
        brdf_value,
        radiance(Ray(hitpoint.position, dir_out),random, depth + 1))
        * cost / pdf;
    return L;
}

};

#endif
