#include "Camera.h"

Camera::Camera(const vec3 &_eye, const vec3 &_center, const vec3 &_up, const float _fovy)
  : eye(_eye)
  , center(_center)
  , up(_up)
  , fovy(_fovy)
{
    
}
