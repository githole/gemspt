#ifndef	_SCENE_H_
#define	_SCENE_H_

#include "constant.h"
#include "sphere.h"
#include "material.h"
#include "hitpoint.h"

namespace gemspt {

#define SCENE_DIFFUSE_ONLY
// #define SCENE_SPECULAR
// #define SCENE_GLASS

class SceneSphere {
private:
    Sphere sphere_;
    const Material *material_;
public:
    SceneSphere(const Sphere &sphere, const Material *material) :
      sphere_(sphere), material_(material) {}

    const Sphere* get_sphere() const {
        return &sphere_;
    }

    const Material* get_material() const {
        return material_;
    }
};

// シーンとの交差判定関数。
inline const SceneSphere* intersect_scene(const Ray &ray, Hitpoint *hitpoint) {
    // レンダリングするシーンデータ。
    // 簡単のため、球のみで構成することにする。
#if defined(SCENE_DIFFUSE_ONLY)
    static const SceneSphere scene[] = {
        SceneSphere(Sphere(100000.0, Vec( 0.0, -100000.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(100000.0, Vec( 0.0,  100004.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(100000.0, Vec(-100003.0,  0.0,       0.0)), new LambertianMaterial(Color(0.7, 0.1, 0.1))),
        SceneSphere(Sphere(100000.0, Vec( 100009.0,  0.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(100000.0, Vec(0.0,        0.0, -100003.0)), new LambertianMaterial(Color(0.1, 0.7, 0.1))),
        SceneSphere(Sphere(100.0,    Vec( 0.0,    103.99,       0.0)), new Lightsource       (Color(8.0, 8.0, 8.0))),
        SceneSphere(Sphere(1.0,      Vec(-2.0,       1.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(1.0,      Vec( 2.0,       1.0,       0.0)), new LambertianMaterial(Color(0.1, 0.1, 0.7))),
    };
#elif defined(SCENE_SPECULAR)
    static const SceneSphere scene[] = {
        SceneSphere(Sphere(100000.0, Vec( 0.0, -100000.0,       0.0)), new PhongMaterial     (Color(0.999, 0.999, 0.999), 100.0)),
        SceneSphere(Sphere(100000.0, Vec( 0.0,  100004.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(100000.0, Vec(-100003.0,  0.0,       0.0)), new LambertianMaterial(Color(0.7, 0.1, 0.1))),
        SceneSphere(Sphere(100000.0, Vec( 100009.0,  0.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(100000.0, Vec(0.0,        0.0, -100003.0)), new LambertianMaterial(Color(0.1, 0.7, 0.1))),
        SceneSphere(Sphere(100.0,    Vec( 0.0,    103.99,       0.0)), new Lightsource       (Color(8.0, 8.0, 8.0))),
        SceneSphere(Sphere(1.0,      Vec(-2.0,       1.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(1.0,      Vec( 2.0,       1.0,       0.0)), new LambertianMaterial(Color(0.1, 0.1, 0.7))),
    };
#elif defined(SCENE_GLASS)
    static const SceneSphere scene[] = {
        SceneSphere(Sphere(100000.0, Vec( 0.0, -100000.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(100000.0, Vec( 0.0,  100004.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(100000.0, Vec(-100003.0,  0.0,       0.0)), new LambertianMaterial(Color(0.7, 0.1, 0.1))),
        SceneSphere(Sphere(100000.0, Vec( 100009.0,  0.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(100000.0, Vec(0.0,        0.0, -100003.0)), new LambertianMaterial(Color(0.1, 0.7, 0.1))),
        SceneSphere(Sphere(100.0,    Vec( 0.0,    103.99,       0.0)), new Lightsource       (Color(8.0, 8.0, 8.0))),
        SceneSphere(Sphere(1.0,      Vec(-2.0,       1.0,       0.0)), new LambertianMaterial(Color(0.7, 0.7, 0.7))),
        SceneSphere(Sphere(1.0,      Vec( 2.0,       1.0,       0.0)), new GlassMaterial     (Color(0.999999, 0.999999, 0.999999), 1.5)),
    };
#endif

    const int n = sizeof(scene) / sizeof(SceneSphere);

    // 初期化
    *hitpoint = Hitpoint();
    const SceneSphere *now_object = NULL;
    
    // 線形探索
    for (int i = 0; i < n; i ++) { 
        Hitpoint tmp_hitpoint;
        if (scene[i].get_sphere()->intersect(ray, &tmp_hitpoint)) {
            if (tmp_hitpoint.distance < hitpoint->distance) {
                *hitpoint = tmp_hitpoint;
                now_object = &scene[i];
            }
        }
    }

    return now_object;
}

};

#endif
