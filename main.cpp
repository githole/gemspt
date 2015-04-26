#include <iostream>
#include "render.h"

int main() {
    std::cout << "gemspt 2015" << std::endl;

    gemspt::render(
        "image.ppm", // 保存ファイル名
        640, 480, // 解像度
        1, // サブピクセルごとのサンプリング数
        4, // サブピクセルの縦横解像度
        8); // スレッド数

    std::cout << "Done." << std::endl;

    return 0;
}