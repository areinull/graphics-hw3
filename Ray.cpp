#include "Ray.h"

Ray::Ray(void)
  : pos(0.0, 0.0, 0.0, 1.0)
  , dir(0.0, 0.0, 1.0, 0.0)
  , t_min(0.0)
  , t_max(99.0)
{
}

Ray& Ray::operator*(const mat4 &M)
{
    pos = pos * M;
//     pos[0] /= pos[3];
//     pos[1] /= pos[3];
//     pos[2] /= pos[3];
    dir = dir * M;
}

bool Ray::operator()(float t, vec4 *v)
{
    if (t < t_min || t > t_max)
        return false;
    *v = pos + dir * t;
    return true;
}
