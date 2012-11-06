#ifndef _LIGHT_H_
#define _LIGHT_H_

#include "Transform.h"

struct Light
{
    enum Light_type {point, directional} ltype;
    vec3 pos;
    vec3 color;
    vec3 attenuation;
};

#endif //_LIGHT_H_