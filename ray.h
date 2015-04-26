#ifndef _RAY_H_
#define _RAY_H_

#include "vec.h"

namespace gemspt {

struct Ray {
    Vec org, dir;
    Ray(const Vec &org, const Vec &dir) : org(org), dir(dir) {}
};

};

#endif
