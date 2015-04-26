#ifndef _PPM_H_
#define _PPM_H_

#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>

namespace gemspt {

inline double clamp(double x){ 
    if (x < 0.0)
        return 0.0;
    if (x > 1.0)
        return 1.0;
    return x;
} 

inline int to_LDR(double x) {
    // ディスプレイのガンマが2.2であることを仮定し、1/2.2乗する。
    // 簡易的なLDR化処理。
    return int(pow(clamp(x), 1/2.2) * 255 + 0.5);
}

void save_ppm_file(const std::string &filename, const Color *image, const int width, const int height) {
    FILE *f = fopen(filename.c_str(), "wb");
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
        fprintf(f,"%d %d %d ", to_LDR(image[i].x), to_LDR(image[i].y), to_LDR(image[i].z));
    fclose(f);
}


};

#endif
