#ifndef _RENDER_H_
#define _RENDER_H_

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "radiance.h"
#include "ppm.h"
#include "random.h"

namespace gemspt {

int render(const char *filename, const int width, const int height, const int num_sample_per_subpixel, const int num_subpixel, const int num_thread) {
#ifdef _OPENMP
    omp_set_num_threads(num_thread);
#endif // _OPENMP
    
    // カメラ位置。
    const Vec camera_position = Vec(7.0, 3.0, 7.0);
    const Vec camera_lookat   = Vec(0.0, 1.0, 0.0);
    const Vec camera_dir      = normalize(camera_lookat - camera_position);
    const Vec camera_up       = Vec(0.0, 1.0, 0.0);

    // ワールド座標系でのイメージセンサーの大きさ。
    const double sensor_width = 30.0 * width / height; // アスペクト比調整。
    const double sensor_height= 30.0;
    // イメージセンサーまでの距離。
    const double sensor_dist  = 45.0;
    // イメージセンサーを張るベクトル。
    const Vec sensor_x_vec = normalize(cross(camera_dir, camera_up)) * sensor_width;
    const Vec sensor_y_vec = normalize(cross(sensor_x_vec, camera_dir)) * sensor_height;
    const Vec sensor_center = camera_position + camera_dir * sensor_dist;

    Color *image = new Color[width * height];
    std::cout << width << "x" << height << " " << num_sample_per_subpixel * (num_subpixel * num_subpixel) << " spp" << std::endl;
    
    for (int y = 0; y < height; ++y) {
        std::cerr << "Rendering (y = " << y << ", " << (100.0 * y / (height - 1)) << " %)          \r";
#pragma omp parallel for schedule(static) // OpenMP
        for (int x = 0; x < width; ++x) {
            Random random(y * width + x + 1);

            const int image_index = (height - y - 1) * width + x;
            // num_subpixel x num_subpixel のスーパーサンプリング。
            for (int sy = 0; sy < num_subpixel; ++sy) {
                for (int sx = 0; sx < num_subpixel; ++sx) {
                    Color accumulated_radiance = Color();
                    // 一つのサブピクセルあたりsamples回サンプリングする。
                    for (int s = 0; s < num_sample_per_subpixel; s ++) {
                        const double rate = (1.0 / num_subpixel);
                        const double r1 = sx * rate + rate / 2.0;
                        const double r2 = sy * rate + rate / 2.0;
                        // イメージセンサー上の位置。
                        const Vec position_on_sensor = 
                            sensor_center + 
                            sensor_x_vec * ((r1 + x) / width - 0.5) +
                            sensor_y_vec * ((r2 + y) / height- 0.5);
                        // レイを飛ばす方向。
                        const Vec dir = normalize(position_on_sensor - camera_position);

                        accumulated_radiance = accumulated_radiance + 
                            radiance(Ray(camera_position, dir), random, 0) 
                            / (double)num_sample_per_subpixel / (double)(num_subpixel * num_subpixel);
                    }
                    image[image_index] = image[image_index] + accumulated_radiance;
                }
            }
        }
    }
    std::cout << std::endl;
    
    // 出力
    save_ppm_file(filename, image, width, height);
    delete[] image;

    return 0;
}


};

#endif
