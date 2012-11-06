#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "Transform.h"

struct Camera
{
    Camera(const vec3 &_eye, const vec3 &_center, const vec3 &_up, const float _fovy);
    
    vec3 eye;    // eye position
    vec3 center; // center look point
    vec3 up;     // up vector
    float fovy;  // field of view [deg]

private:
    Camera(void);
};

#endif // _CAMERA_H_