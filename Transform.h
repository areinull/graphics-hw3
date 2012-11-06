#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using glm::mat3;
using glm::mat4; 
using glm::vec3; 
using glm::vec4; 
const float pi = 3.14159265 ; // For portability across platforms

namespace Transform  
{
    void left(float degrees, vec3& eye, vec3& up);
    void up(float degrees, vec3& eye, vec3& up);
    mat3 rotate(const float degrees, const vec3& axis);
    mat4 scale(const float &sx, const float &sy, const float &sz);
    mat4 translate(const float &tx, const float &ty, const float &tz);
    vec3 upvector(const vec3 &up, const vec3 &zvec);
};

#endif // _TRANSFORM_H_
