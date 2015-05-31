#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include <assert.h>

#include "random.h"
#include "vec.h"
#include "constant.h"
#include "sampling.h"

namespace gemspt {

typedef Vec Color;

// マテリアルインターフェース
class Material {
protected:
    Color emission_;
    Color reflectance_;
public:
    Material(const Color &emission, const Color &reflectance) : emission_(emission), reflectance_(reflectance) {}
    virtual Color emission() const {
        return emission_;
    }
    virtual Color reflectance() const {
        return reflectance_;
    }

    // in, outはカメラ側から光を逆方向に追跡したときの入出方向とする。
    // 以下、in = -omega, out = omega'となる。
    virtual Color eval(const Vec &in, const Vec &normal, const Vec &out) const = 0; // BRDFとして評価した時の値。
    virtual Vec sample(Random &random, const Vec &in, const Vec &normal, double *pdf, Color *brdf_value) const = 0; // 次の反射方向をサンプリング。
};

// Lambertian BRDF
// 所謂完全拡散面
class LambertianMaterialSimple : public Material {
private:
public:
    LambertianMaterialSimple(const Color &reflectance) : Material(Color(), reflectance) {};

    // Lambertian BRDFはρ/πになる（ρは反射率）
    virtual Color eval(const Vec &in, const Vec &normal, const Vec &out) const {
        return reflectance_ / kPI;
    }
    
    // 単純に半球一様サンプリングする。
    virtual Vec sample(Random &random, const Vec &in, const Vec &normal, double *pdf, Color *brdf_value) const {
        Vec binormal, tangent, now_normal = normal;

        createOrthoNormalBasis(now_normal, &tangent, &binormal);
        const Vec dir = Sampling::uniformHemisphereSurface(random, now_normal, tangent, binormal);

        // pdf: 1/(2 * pi)
        if (pdf != NULL) {
            *pdf = 1.0 / (2.0 * kPI);
        }
        if (brdf_value != NULL) {
            *brdf_value = eval(in, normal, dir);
        }
        return dir;
    }
};

// Lambertian BRDF
// 所謂完全拡散面
// インポータンスサンプリング版。
class LambertianMaterial : public Material {
private:
public:
    LambertianMaterial(const Color &reflectance) : Material(Color(), reflectance) {};

    // Lambertian BRDFはρ/πになる（ρは反射率）
    virtual Color eval(const Vec &in, const Vec &normal, const Vec &out) const {
        return reflectance_ / kPI;
    }

    // pdfとしてcosΘ/piを使用してインポータンスサンプリングする。
    virtual Vec sample(Random &random, const Vec &in, const Vec &normal, double *pdf, Color *brdf_value) const {
        Vec binormal, tangent, now_normal = normal;

        createOrthoNormalBasis(now_normal, &tangent, &binormal);
        const Vec dir = Sampling::cosineWeightedHemisphereSurface(random, now_normal, tangent, binormal);

        // pdf: cosΘ/pi
        if (pdf != NULL) {
            *pdf = dot(normal, dir) / kPI;
        }
        if (brdf_value != NULL) {
            *brdf_value = eval(in, normal, dir);
        }
        return dir;
    }
};

// 正規化Phong BRDF
class PhongMaterial : public Material {
private:
    double n_;
public:
    PhongMaterial(const Color &reflectance, const double n) : Material(Color(), reflectance), n_(n) {}

    // 正規化Phong
    virtual Color eval(const Vec &in, const Vec &normal, const Vec &out) const {
        if (dot(normal, out) < 0) {
            // 次のレイの方向が地面より下の方向だったら0。
            return Color();
        }

        const Vec reflection_dir = reflect(in, normal);
        double cosa = dot(reflection_dir, out);
        if (cosa < 0)
            cosa = 0.0;
        return reflectance_ * (n_ + 2.0) / (2.0 * kPI) * pow(cosa, n_);
    }

    // BRDF形状をpdfとして使ってインポータンスサンプリングする。
    virtual Vec sample(Random &random, const Vec &in, const Vec &normal, double *pdf, Color *brdf_value) const {
        Vec dir;
        const Vec reflection_dir = reflect(in, normal);
        Vec binormal, tangent;
        createOrthoNormalBasis(reflection_dir, &tangent, &binormal);

        const double u1 = random.next01();
        const double u2 = random.next01();
        
        const double phi = u1 * 2.0 * kPI;
        const double theta = acos(pow(u2, 1 / (n_ + 1)));

        dir = tangent * sin(theta) * cos(phi) + reflection_dir * cos(theta) + binormal *sin(theta) * sin(phi);
        
        if (pdf != NULL) {
            double cosa = dot(reflection_dir, dir);
            if (cosa < 0)
                cosa = 0.0;
            *pdf = (n_ + 1.0) / (2.0 * kPI) * pow(cosa, n_);
        }
        if (brdf_value != NULL) {
            *brdf_value = eval(in, normal, dir);
        }

        return dir;
    }
};

// 理想的なガラス面。
#define DELTA 1.0
class GlassMaterial : public Material {
private:
    double ior_;
public:
    GlassMaterial(const Color &reflectance, const double ior) : Material(Color(), reflectance), ior_(ior) {}
    
    // 理想的なガラス面におけるBRDFはディラックのδ関数を使ってδ/cosΘとなる。
    // δ関数を表現することは出来ないが、モンテカルロ積分においてはpdfにもδ関数が現れるため分母と分子で打ち消し合う。
    // そこでcosΘと反射率だけ入れておく。
    // in, outはカメラ側から追跡したときの入出方向なので、光の入射方向はoutになるため、cosθは法線とoutの内積になる。
    // 反射率や透過率（Fr,Ft）は含まれていないことに注意！
    virtual Color eval(const Vec &in, const Vec &normal, const Vec &out) const {
        return reflectance_ * DELTA / dot(normal, out);
    }

    virtual Vec sample(Random &random, const Vec &in, const Vec &normal, double *pdf, Color *brdf_value) const {
        const Vec now_normal = dot(normal, in) < 0.0 ? normal: -normal; // 交差位置の法線（物体からのレイの入出を考慮。
        const bool into = dot(normal, now_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか。
        const double n1 = 1.0; // 真空の屈折率
        const double n2 = ior_; // オブジェクトの屈折率
        const double n = into ? n1 / n2 : n2 / n1;

        // Snellの法則を用いて屈折方向を計算する。
        const double dir_dot_normal = dot(in, now_normal);
        const double cos2t_2 = 1.0 - n * n * (1.0 - dir_dot_normal * dir_dot_normal);
        
        // 全反射
        const Vec reflection_dir = reflect(in, now_normal);
        if (cos2t_2 < 0.0) {
            if (pdf != NULL) {
                // pdfはディラックのδ関数なので実数値にはならないが、将来的にモンテカルロ積分において、
                // 分母と分子の両方にδが表れるため結局打ち消し合うため、1でよい。あくまでδであること忘れないためにDELTAを入れておくが、実態は1。
                *pdf = DELTA;
            }
            if (brdf_value != NULL) {
                *brdf_value = eval(in, normal, reflection_dir);
            }
            return reflection_dir;
        }

        // 屈折の方向
        const Vec refraction_dir = in * n - now_normal * (dir_dot_normal * n + sqrt(cos2t_2));
        
        // Fresnelの式
        const double cost_1 = dot(-in, now_normal);
        const double cost_2 = sqrt(cos2t_2);
        const double r_parallel = (n * cost_1 - cost_2) / (n * cost_1 + cost_2);
        const double r_perpendicular = (cost_1 - n * cost_2) / (cost_1 + n * cost_2);
        const double Fr = 0.5 * (r_parallel * r_parallel + r_perpendicular * r_perpendicular);

        const double factor = pow(into ? n1 / n2 : n2 / n1, 2.0); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
        const double Ft = (1.0 - Fr) * factor; // 屈折方向の割合。
        
        // ロシアンルーレットで屈折か反射かを決定する。
        // ロシアンルーレットの確率は反射率ということしておく。
        const double probability  = Fr;
        if (random.next01() < probability) { // 反射
            if (pdf != NULL) {
                // pdfはディラックのδ関数なので実数値にはならないが、将来的にモンテカルロ積分において、
                // 分母と分子の両方にδが表れるため結局打ち消し合うため、1でよい。あくまでδであること忘れないためにDELTAを入れておくが、実態は1。
                *pdf = DELTA * probability;
            }
            if (brdf_value != NULL) {
                *brdf_value = Fr * eval(in, normal, reflection_dir);
            }
            return reflection_dir;
        } else { // 屈折
            if (pdf != NULL) {
                *pdf = DELTA * (1.0f - probability);
            }
            if (brdf_value != NULL) {
                *brdf_value = Ft * eval(in, normal, reflection_dir);
            }
            return refraction_dir;
        }
    }
};
#undef DELTA

// 光源としてふるまうマテリアル
class Lightsource : public Material {
private:
public:
    Lightsource(const Color &emission) : Material(emission, Color()) {}

    virtual Color eval(const Vec &in, const Vec &normal, const Vec &out) const {
        assert(false);
        return Color();
    }

    virtual double eval_pdf(const Vec &in, const Vec &normal, const Vec &out) const {
        assert(false);
        return -1;
    }

    virtual Vec sample(Random &random, const Vec &in, const Vec &normal, double *pdf, Color *brdf_value) const {
        assert(false);
        return Color();
    }
};

};

#endif
