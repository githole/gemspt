#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "vec.h"
#include "constant.h"

namespace gemspt {

struct Hitpoint {
    double distance;
    Vec normal;
    Vec position;

    Hitpoint() : distance(kINF), normal(), position() {}
};

};

#endif
