#ifndef _BOUNDINGVOLUME_H_
#define _BOUNDINGVOLUME_H_

#include "Transform.h"
#include <memory>
#include <deque>

class BoundingVolume;
typedef std::shared_ptr<BoundingVolume> BVp;

// axis aligned bounding volume box
class BoundingVolume
{
    friend class Scene;

    vec3 lower; // lower corner of bv
    vec3 upper; // upper corner of bv
    BVp child1, child2;
    void *obj;
    
public:
    BoundingVolume(): lower(HUGE_VALF), upper(HUGE_VALF) {}
    
    void expand(BVp bv);
    void addChild(BVp bv);
    
    static void sort(std::deque<BVp> &bvl, int start, int end, int axis);
};

#endif // _BOUNDINGVOLUME_H_
