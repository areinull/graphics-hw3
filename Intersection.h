#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "Transform.h"

struct Intersection
{
    bool hasIntersection;
    vec4 n;
    float t;
    vec4 pos;
    vec4 sourcedir;
    const void *obj; // hit object
};

#endif // _INTERSECTION_H_