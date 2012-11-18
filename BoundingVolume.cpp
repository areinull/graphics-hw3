#include "BoundingVolume.h"
#include <algorithm>

void BoundingVolume::expand(BVp bv)
{
    if (lower[0] == HUGE_VALF)
    {
        lower = bv->lower;
        upper = bv->upper;
        return;
    }
    for (int i = 0; i < 3; i++)
    {
        if (bv->lower[i] < lower[i])
            lower[i] = bv->lower[i];
        if (bv->upper[i] > upper[i])
            upper[i] = bv->upper[i];
    }
}

void BoundingVolume::addChild(BVp bv)
{
    if (!child1)
        child1 = bv;
    else if (!child2)
        child2 = bv;
    else
        throw("Requested to add child while have 2 already");
}

void BoundingVolume::sort(std::deque< BVp >& bvl, int start, int end, int axis)
{
    std::deque< BVp >::iterator s = bvl.begin();
    std::deque< BVp >::iterator e = bvl.begin();
    std::advance(s, start);
    std::advance(e, end);
    std::sort(s, e, [axis] (BVp lhs, BVp rhs) -> bool
        {
            return lhs->lower[axis] < rhs->lower[axis];
        });
}
