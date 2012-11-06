#ifndef _RAY_H_
#define _RAY_H_

#include "Transform.h"

struct Ray
{
    Ray(void);
    Ray& operator*(const mat4 &M);
    bool operator()(float t, vec4 *v);
    
    vec4 pos;    // position
    vec4 dir;    // direction
    float t_min, t_max;
};


#endif // _RAY_H_