#include "Transform.h"

mat3 Transform::rotate(const float degrees, const vec3& axis) 
{
    vec3 normAxis = glm::normalize(axis);
    mat3 R;
    float rad = degrees * pi / 180.f;
  
    R = cos(rad)*mat3(1, 0, 0,
                      0, 1, 0,
                      0, 0, 1);
    mat3 tmp;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            tmp[i][j] = normAxis[i]*normAxis[j];
    R += (1.-cos(rad)) * tmp;
    R += sin(rad) * mat3(0, -normAxis[2], normAxis[1],
                         normAxis[2], 0, -normAxis[0],
                         -normAxis[1], normAxis[0], 0);
    return R;
}

void Transform::left(float degrees, vec3& eye, vec3& up) 
{ 
    mat3 R = rotate(degrees, up);
    eye = eye * R;
}

void Transform::up(float degrees, vec3& eye, vec3& up) 
{
    vec3 aux = glm::normalize(glm::cross(eye, up));
    mat3 R = rotate(degrees, aux);
    eye = eye * R;
    up = up * R;
}

mat4 Transform::scale(const float &sx, const float &sy, const float &sz) 
{
    return mat4(sx, 0, 0, 0,
                0, sy, 0, 0,
                0, 0, sz, 0,
                0, 0, 0,  1);
}

mat4 Transform::translate(const float &tx, const float &ty, const float &tz) 
{
    return mat4(1, 0, 0, tx,
                0, 1, 0, ty,
                0, 0, 1, tz,
                0, 0, 0, 1);
}

vec3 Transform::upvector(const vec3 &up, const vec3 & zvec) 
{
    vec3 x = glm::cross(up,zvec); 
    vec3 y = glm::cross(zvec,x); 
    vec3 ret = glm::normalize(y); 
    return ret; 
}
