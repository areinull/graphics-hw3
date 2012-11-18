#include "Scene.h"
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <glm/gtx/color_cast.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <boost/thread.hpp>
#include <list>

Scene::Scene(string fname)
  : w(160)
  , h(120)
  , aspect(4.f/3.f)
  , maxdepth(5)
  , fb(0)
  , cam(vec3(0.0, 0.0, 5.0), vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), 90.0)
  , saveFile("screenshot.png")
  , maxverts(0)
  , usedverts(0)
  , vertex(0)
  , maxvertnorms(0)
  , usedvertnorms(0)
  , vertnorm(0)
{
    readfile(fname);
    createBV();
}

Scene::~Scene(void)
{
    delete [] fb;
    delete [] vertex;
    delete [] vertnorm;
    for(vector<Object*>::iterator it=objs.begin(); it != objs.end(); ++it)
        delete *it;
}

void Scene::saveScreenshot(void)
{
    FIBITMAP *img = FreeImage_ConvertFromRawBits(fb, w, h, w*3, 24, 0xFF0000, 0x00FF00, 0x0000FF, true);
    cout << "Saving screenshot: " << saveFile << "\n";
    FreeImage_Save(FIF_PNG, img, saveFile.c_str(), 0);
}

void Scene::worker(int x)
{
    for (int i = 0; i < h; i++)
    {
        for (int j = x; j < w; j+=8)
        {
            Ray r = getRay(j, i);
            Intersection hit = intersect(r);
            vec3 color = findColor(hit, 0);
            fb[3*w*i + 3*j + 2] = glm::u8channel_cast(color.r);   // red
            fb[3*w*i + 3*j + 1] = glm::u8channel_cast(color.g);   // green
            fb[3*w*i + 3*j + 0] = glm::u8channel_cast(color.b);   // blue
        }
        if (!x)
            cout << '\r' << 100*i/h << "%" << flush;
    }
}

static void call(Scene *s, int x)
{
    s->worker(x);
}

void Scene::render(void)
{
    using namespace boost;
    fb = new BYTE[3*w*h];
    aspect = (float)w / (float)h;
    list<thread*> threads;
    for (int i=0; i<8; i++)
    {
        threads.push_back(new thread(call, this, i));
    }
    for (list<thread*>::iterator it=threads.begin(); it!=threads.end(); ++it)
    {
        (*it)->join();
    }
    for (list<thread*>::iterator it=threads.begin(); it!=threads.end(); ++it)
    {
        delete *it;
    }
    threads.clear();
    cout << '\r' << "Done" << endl;
}

vec3 Scene::findColor(const Intersection& hit, int depth)
{
    vec3 color(0.0,0.0,0.0);
    if (hit.hasIntersection)
    {
        const Object &obj = *(const Object*)hit.obj;
        color += obj.ambient + obj.emission;
        for (vector<Light>::const_iterator l = lights.begin(); l != lights.end(); ++l)
        {
            Ray ray;
            ray.pos = hit.pos + 0.0002*hit.n;
            ray.dir = glm::normalize(vec4(l->pos,1.0) - ray.pos);
            Intersection block = intersect(ray);
            float dist = glm::distance(vec3(hit.pos[0],hit.pos[1],hit.pos[2]), l->pos);
            if (block.hasIntersection)
            {
                if (l->ltype == Light::directional || (l->ltype == Light::point && dist > block.t))
                    continue;
            }
            vec3 L = l->color;
            if (l->ltype == Light::point)
                L /= l->attenuation[0]+l->attenuation[1]*dist+l->attenuation[2]*dist*dist;
            vec4 direction;
            if (l->ltype == Light::directional)
                direction = vec4(glm::normalize(l->pos),0.0);
            else
                direction = glm::normalize(vec4(l->pos,1.0) - hit.pos);
            vec4 half = glm::normalize(direction + hit.sourcedir);
            vec3 lambert = obj.diffuse * max(glm::dot(hit.n,direction),0.0f);
            vec3 phong = obj.specular * pow(max(glm::dot(hit.n,half),0.0f), obj.shininess);
            color += L * (lambert + phong);
        }
        if (depth < maxdepth)
        {
            Ray ray;
            ray.pos = hit.pos + 0.0002*hit.n;
            ray.dir = glm::normalize(2.*glm::dot(hit.sourcedir, hit.n)*hit.n - hit.sourcedir); // mirror direction
            Intersection nextHit = intersect(ray);
            color += obj.specular * findColor(nextHit, depth+1);
        }
    }
    return color;
}

bool Scene::intersectBV(Ray& r, BVp bv)
{
    /* 
     Fast Ray-Box Intersection
     by Andrew Woo
     from "Graphics Gems", Academic Press, 1990
    */
    const int RIGHT = 0;
    const int LEFT = 1;
    const int MIDDLE = 2;
    vec3 origin(r.pos[0]/r.pos[3], r.pos[1]/r.pos[3], r.pos[2]/r.pos[3]);
    vec3 coord; // hit point
    bool inside = true;
    int quadrant[3];
    register int i;
    int whichPlane;
    vec3 maxT;
    vec3 candidatePlane;
    
    /* Find candidate planes; this loop can be avoided if
     * rays cast all from the eye(assume perpsective view) */
    for (i = 0; i < 3; i++)
    {
        if(origin[i] < bv->lower[i])
        {
            quadrant[i] = LEFT;
            candidatePlane[i] = bv->lower[i];
            inside = false;
        }
        else if (origin[i] > bv->upper[i])
        {
            quadrant[i] = RIGHT;
            candidatePlane[i] = bv->upper[i];
            inside = false;
        }
        else
        {
            quadrant[i] = MIDDLE;
        }
    }
    
    /* Ray origin inside bounding box */
    if (inside)
    {
        return true;
    }
    
    
    /* Calculate T distances to candidate planes */
    for (i = 0; i < 3; i++)
    {
        if (quadrant[i] != MIDDLE && r.dir[i] != 0.)
            maxT[i] = (candidatePlane[i] - origin[i]) / r.dir[i];
        else
            maxT[i] = -1.;
    }
        
    /* Get largest of the maxT's for final choice of intersection */
    whichPlane = 0;
    for (i = 1; i < 3; i++)
        if (maxT[whichPlane] < maxT[i])
            whichPlane = i;
        
    /* Check final candidate actually inside box */
    if (maxT[whichPlane] < 0.)
        return false;
    for (i = 0; i < 3; i++)
    {
        if (whichPlane != i)
        {
            coord[i] = origin[i] + maxT[whichPlane] * r.dir[i];
            if (coord[i] < bv->lower[i] || coord[i] > bv->upper[i])
                return false;
        }
        else
        {
            coord[i] = candidatePlane[i];
        }
    }
    return true;      /* ray hits box */
}


Intersection Scene::intersect(Ray& r)
{
    float mindist = HUGE_VALF;
    Intersection i;
    i.hasIntersection = false;
    
//     vector<Object*> cand;
    vector<BVp> stack;
    stack.push_back(rootBV);
    while (!stack.empty())
    {
        BVp cur = stack.back();
        stack.pop_back();
        
        // check for intersection with BV
        if (!intersectBV(r, cur))
            continue;
        
        // if it's the leaf node then calculate intersection
        if (!cur->child1 && !cur->child2)
        {
            Object *obj = (Object*)cur->obj;
            switch (obj->type)
            {
                case tri:
                {
                    vec4 norm = vec4(obj->vert[3], 0.0);
                    if (glm::dot(r.dir, norm) == 0)
                        break;
                    float t = (glm::dot(vec4(obj->vert[0],1.0),norm)-glm::dot(r.pos, norm))/glm::dot(r.dir, norm);
                    if (t < r.t_min || t > r.t_max || t > mindist)
                        break;
                    vec4 cross4;
                    if (!r(t, &cross4))
                        break;
                    vec3 cross3(cross4[0]/cross4[3],cross4[1]/cross4[3],cross4[2]/cross4[3]);
                    if (glm::dot(glm::cross(obj->vert[1]-obj->vert[0],cross3-obj->vert[0]),obj->vert[3]) >= 0 &&
                        glm::dot(glm::cross(obj->vert[2]-obj->vert[1],cross3-obj->vert[1]),obj->vert[3]) >= 0 &&
                        glm::dot(glm::cross(obj->vert[0]-obj->vert[2],cross3-obj->vert[2]),obj->vert[3]) >= 0)
                    {
                        i.hasIntersection = true;
                        i.n = norm;
                        i.t = t;
                        i.pos = vec4(cross3, 1.0);
                        i.sourcedir = -r.dir;
                        i.obj = obj;
                        mindist = t;
                    }
                    break;
                }
                case trinormal:
                {
                    vec3 norm3 = glm::normalize(glm::cross(obj->vert_norm[1].v-obj->vert_norm[0].v,obj->vert_norm[2].v-obj->vert_norm[0].v));
                    vec4 norm = vec4(norm3, 0.0);
                    if (glm::dot(r.dir, norm) == 0)
                        break;
                    float t = (glm::dot(vec4(obj->vert_norm[0].v,1.0),norm)-glm::dot(r.pos, norm))/glm::dot(r.dir, norm);
                    if (t < r.t_min || t > r.t_max || t > mindist)
                        break;
                    vec4 cross4;
                    if (!r(t, &cross4))
                        break;
                    vec3 cross3(cross4[0]/cross4[3],cross4[1]/cross4[3],cross4[2]/cross4[3]);
                    if (glm::dot(glm::cross(obj->vert_norm[1].v-obj->vert_norm[0].v,cross3-obj->vert_norm[0].v),norm3) >= 0 &&
                        glm::dot(glm::cross(obj->vert_norm[2].v-obj->vert_norm[1].v,cross3-obj->vert_norm[1].v),norm3) >= 0 &&
                        glm::dot(glm::cross(obj->vert_norm[0].v-obj->vert_norm[2].v,cross3-obj->vert_norm[2].v),norm3) >= 0)
                    {
                        i.hasIntersection = true;
                        i.t = t;
                        i.pos = vec4(cross3, 1.0);
                        i.sourcedir = -r.dir;
                        i.obj = obj;
                        mindist = t;
                        // interpolate normal at intersection
                        // calculate distances^2 from intersection to choose farthest
                        float IA = glm::dot(obj->vert_norm[0].v-cross3,obj->vert_norm[0].v-cross3);
                        float IB = glm::dot(obj->vert_norm[1].v-cross3,obj->vert_norm[1].v-cross3);
                        float IC = glm::dot(obj->vert_norm[2].v-cross3,obj->vert_norm[2].v-cross3);
                        vec3 A, An; // farthest
                        vec3 B, Bn, C, Cn; // another two
                        if (IA >= IB && IA >= IC)
                        {
                            A = obj->vert_norm[0].v;
                            B = obj->vert_norm[1].v;
                            C = obj->vert_norm[2].v;
                            An = obj->vert_norm[0].n;
                            Bn = obj->vert_norm[1].n;
                            Cn = obj->vert_norm[2].n;
                        }
                        else if (IB >= IA && IB >= IC)
                        {
                            A = obj->vert_norm[1].v;
                            B = obj->vert_norm[2].v;
                            C = obj->vert_norm[0].v;
                            An = obj->vert_norm[1].n;
                            Bn = obj->vert_norm[2].n;
                            Cn = obj->vert_norm[0].n;
                        }
                        else
                        {
                            A = obj->vert_norm[2].v;
                            B = obj->vert_norm[0].v;
                            C = obj->vert_norm[1].v;
                            An = obj->vert_norm[2].n;
                            Bn = obj->vert_norm[0].n;
                            Cn = obj->vert_norm[1].n;
                        }
                        // find Q - intersection of AI and BC
                        Ray AI;
                        AI.pos = vec4(A,1.0);
                        AI.dir = vec4(cross3-A,0.0);
                        // form normal to plane containing BC
                        vec4 norm_aux(glm::normalize(glm::cross(glm::cross(C-B, cross3-A), C-B)),0.0);
                        float AQ_d = (glm::dot(vec4(B,1.0),norm_aux)-glm::dot(AI.pos, norm_aux))/glm::dot(AI.dir, norm_aux);
                        // find Q
                        vec4 Q4;
                        AI(AQ_d, &Q4);
                        vec3 Q(Q4[0]/Q4[3],Q4[1]/Q4[3],Q4[2]/Q4[3]);
                        // interpolate normal at Q
                        float BQ_d = sqrt(glm::dot(B-Q,B-Q));
                        float CQ_d = sqrt(glm::dot(C-Q,C-Q));
                        float BC_d = sqrt(glm::dot(B-C,B-C));
                        vec3 Q_norm = (Bn*CQ_d + Cn*BQ_d)/BC_d;
                        // interpolate normal at I
                        float AI_d = sqrt(glm::dot(A-cross3,A-cross3));
                        float QI_d = sqrt(glm::dot(Q-cross3,Q-cross3));
                        i.n = vec4(glm::normalize((Q_norm*AI_d + An*QI_d)/AQ_d), 0.0);
                    }
                    break;
                }
                case sphere:
                {
                    Ray r1 = r;
                    r1 * glm::inverse(obj->transform);
                    vec4 pos = obj->pos_r;
                    float rad = pos[3];
                    pos[3] = 1.0;
                    float D = pow(glm::dot(r1.dir, r1.pos-pos),2) - glm::dot(r1.dir,r1.dir)*(glm::dot(r1.pos-pos,r1.pos-pos)-rad*rad);
                    if (D < 0.)
                        break;
                    float t, t1, t2;
                    t1 = (-glm::dot(r1.dir, r1.pos-pos)+sqrt(D))/glm::dot(r1.dir,r1.dir);
                    t2 = (-glm::dot(r1.dir, r1.pos-pos)-sqrt(D))/glm::dot(r1.dir,r1.dir);
                    if (t1 < 0)
                        t = t2;
                    else if (t2 < 0)
                        t = t1;
                    else if (t1 > t2)
                        t = t2;
                    else
                        t = t1;
                    if (t < r1.t_min || t > r1.t_max || t > mindist)
                        break;
                    vec4 cross4;
                    if (!r1(t, &cross4))
                        break;
                    i.hasIntersection = true;
                    // some problems with translation component
                    i.n = glm::normalize(cross4 - pos) * glm::inverseTranspose(obj->transform);
                    i.n[3] = 0.0;
                    i.n = glm::normalize(i.n);
                    i.pos = cross4 * obj->transform;
                    i.pos[0] /= i.pos[3];
                    i.pos[1] /= i.pos[3];
                    i.pos[2] /= i.pos[3];
                    i.pos[3] = 1;
                    i.t = glm::distance(vec3(r.pos[0],r.pos[1],r.pos[2]), vec3(i.pos[0],i.pos[1],i.pos[2]));
                    i.sourcedir = -r.dir;
                    i.obj = obj;
                    mindist = i.t;
                    break;
                }
            }   // switch
        }
        else
        {
            if (cur->child2)
                stack.push_back(cur->child2);
            if (cur->child1)
                stack.push_back(cur->child1);
        }
    }
    /*
    for (vector<Object*>::const_iterator it = objs.begin(); it != objs.end(); ++it)
    {
        switch ((*it)->type)
        {
            case tri:
            {
                vec4 norm = vec4((*it)->vert[3], 0.0);
                if (glm::dot(r.dir, norm) == 0)
                    break;
                float t = (glm::dot(vec4((*it)->vert[0],1.0),norm)-glm::dot(r.pos, norm))/glm::dot(r.dir, norm);
                if (t < r.t_min || t > r.t_max || t > mindist)
                    break;
                vec4 cross4;
                if (!r(t, &cross4))
                    break;
                vec3 cross3(cross4[0]/cross4[3],cross4[1]/cross4[3],cross4[2]/cross4[3]);
                if (glm::dot(glm::cross((*it)->vert[1]-(*it)->vert[0],cross3-(*it)->vert[0]),(*it)->vert[3]) >= 0 &&
                    glm::dot(glm::cross((*it)->vert[2]-(*it)->vert[1],cross3-(*it)->vert[1]),(*it)->vert[3]) >= 0 &&
                    glm::dot(glm::cross((*it)->vert[0]-(*it)->vert[2],cross3-(*it)->vert[2]),(*it)->vert[3]) >= 0)
                {
                    i.hasIntersection = true;
                    i.n = norm;
                    i.t = t;
                    i.pos = vec4(cross3, 1.0);
                    i.sourcedir = -r.dir;
                    i.obj = *it;
                    mindist = t;
                }
                break;
            }
            case trinormal:
            {
                vec3 norm3 = glm::normalize(glm::cross((*it)->vert_norm[1].v-(*it)->vert_norm[0].v,(*it)->vert_norm[2].v-(*it)->vert_norm[0].v));
                vec4 norm = vec4(norm3, 0.0);
                if (glm::dot(r.dir, norm) == 0)
                    break;
                float t = (glm::dot(vec4((*it)->vert_norm[0].v,1.0),norm)-glm::dot(r.pos, norm))/glm::dot(r.dir, norm);
                if (t < r.t_min || t > r.t_max || t > mindist)
                    break;
                vec4 cross4;
                if (!r(t, &cross4))
                    break;
                vec3 cross3(cross4[0]/cross4[3],cross4[1]/cross4[3],cross4[2]/cross4[3]);
                if (glm::dot(glm::cross((*it)->vert_norm[1].v-(*it)->vert_norm[0].v,cross3-(*it)->vert_norm[0].v),norm3) >= 0 &&
                    glm::dot(glm::cross((*it)->vert_norm[2].v-(*it)->vert_norm[1].v,cross3-(*it)->vert_norm[1].v),norm3) >= 0 &&
                    glm::dot(glm::cross((*it)->vert_norm[0].v-(*it)->vert_norm[2].v,cross3-(*it)->vert_norm[2].v),norm3) >= 0)
                {
                    i.hasIntersection = true;
                    i.t = t;
                    i.pos = vec4(cross3, 1.0);
                    i.sourcedir = -r.dir;
                    i.obj = *it;
                    mindist = t;
                    // interpolate normal at intersection
                    // calculate distances^2 from intersection to choose farthest
                    float IA = glm::dot((*it)->vert_norm[0].v-cross3,(*it)->vert_norm[0].v-cross3);
                    float IB = glm::dot((*it)->vert_norm[1].v-cross3,(*it)->vert_norm[1].v-cross3);
                    float IC = glm::dot((*it)->vert_norm[2].v-cross3,(*it)->vert_norm[2].v-cross3);
                    vec3 A, An; // farthest
                    vec3 B, Bn, C, Cn; // another two
                    if (IA >= IB && IA >= IC)
                    {
                        A = (*it)->vert_norm[0].v;
                        B = (*it)->vert_norm[1].v;
                        C = (*it)->vert_norm[2].v;
                        An = (*it)->vert_norm[0].n;
                        Bn = (*it)->vert_norm[1].n;
                        Cn = (*it)->vert_norm[2].n;
                    }
                    else if (IB >= IA && IB >= IC)
                    {
                        A = (*it)->vert_norm[1].v;
                        B = (*it)->vert_norm[2].v;
                        C = (*it)->vert_norm[0].v;
                        An = (*it)->vert_norm[1].n;
                        Bn = (*it)->vert_norm[2].n;
                        Cn = (*it)->vert_norm[0].n;
                    }
                    else
                    {
                        A = (*it)->vert_norm[2].v;
                        B = (*it)->vert_norm[0].v;
                        C = (*it)->vert_norm[1].v;
                        An = (*it)->vert_norm[2].n;
                        Bn = (*it)->vert_norm[0].n;
                        Cn = (*it)->vert_norm[1].n;
                    }
                    // find Q - intersection of AI and BC
                    Ray AI;
                    AI.pos = vec4(A,1.0);
                    AI.dir = vec4(cross3-A,0.0);
                    // form normal to plane containing BC
                    vec4 norm_aux(glm::normalize(glm::cross(glm::cross(C-B, cross3-A), C-B)),0.0);
                    float AQ_d = (glm::dot(vec4(B,1.0),norm_aux)-glm::dot(AI.pos, norm_aux))/glm::dot(AI.dir, norm_aux);
                    // find Q
                    vec4 Q4;
                    AI(AQ_d, &Q4);
                    vec3 Q(Q4[0]/Q4[3],Q4[1]/Q4[3],Q4[2]/Q4[3]);
                    // interpolate normal at Q
                    float BQ_d = sqrt(glm::dot(B-Q,B-Q));
                    float CQ_d = sqrt(glm::dot(C-Q,C-Q));
                    float BC_d = sqrt(glm::dot(B-C,B-C));
                    vec3 Q_norm = (Bn*CQ_d + Cn*BQ_d)/BC_d;
                    // interpolate normal at I
                    float AI_d = sqrt(glm::dot(A-cross3,A-cross3));
                    float QI_d = sqrt(glm::dot(Q-cross3,Q-cross3));
                    i.n = vec4(glm::normalize((Q_norm*AI_d + An*QI_d)/AQ_d), 0.0);
                }
                break;
            }
            case sphere:
            {
                Ray r1 = r;
                r1 * glm::inverse((*it)->transform);
                vec4 pos = (*it)->pos_r;
                float rad = pos[3];
                pos[3] = 1.0;
                float D = pow(glm::dot(r1.dir, r1.pos-pos),2) - glm::dot(r1.dir,r1.dir)*(glm::dot(r1.pos-pos,r1.pos-pos)-rad*rad);
                if (D < 0.)
                    break;
                float t, t1, t2;
                t1 = (-glm::dot(r1.dir, r1.pos-pos)+sqrt(D))/glm::dot(r1.dir,r1.dir);
                t2 = (-glm::dot(r1.dir, r1.pos-pos)-sqrt(D))/glm::dot(r1.dir,r1.dir);
                if (t1 < 0)
                    t = t2;
                else if (t2 < 0)
                    t = t1;
                else if (t1 > t2)
                    t = t2;
                else
                    t = t1;
                if (t < r1.t_min || t > r1.t_max || t > mindist)
                    break;
                vec4 cross4;
                if (!r1(t, &cross4))
                    break;
                i.hasIntersection = true;
                // some problems with translation component
                i.n = glm::normalize(cross4 - pos) * glm::inverseTranspose((*it)->transform);
                i.n[3] = 0.0;
                i.n = glm::normalize(i.n);
                i.pos = cross4 * (*it)->transform;
                i.pos[0] /= i.pos[3];
                i.pos[1] /= i.pos[3];
                i.pos[2] /= i.pos[3];
                i.pos[3] = 1;
                i.t = glm::distance(vec3(r.pos[0],r.pos[1],r.pos[2]), vec3(i.pos[0],i.pos[1],i.pos[2]));
                i.sourcedir = -r.dir;
                i.obj = *it;
                mindist = i.t;
                break;
            }
        }
    }
    */
    return i;
}

Ray Scene::getRay(int i, int j)
{
    Ray r;
    r.pos = vec4(cam.eye, 1.0);
    vec3 w = glm::normalize(cam.eye - cam.center);
    vec3 u = glm::cross(cam.up, w);
    vec3 v = glm::cross(w, u);
    float alpha = tanf(cam.fovy*pi/360.f)*aspect*(i + 0.5f - this->w/2.f)/(this->w/2.f);
    float beta = tanf(cam.fovy*pi/360.f)*(h/2.f - j - 0.5f)/(h/2.f);
    vec3 dir3 = alpha*u + beta*v - w;
    r.dir = vec4(glm::normalize(dir3), 0.0);
    return r;
}

// The function below applies the appropriate transform to a 4-vector
void Scene::matransform(stack<mat4> &transfstack, float* values)
{
    mat4 transform = transfstack.top();
    vec4 valvec = vec4(values[0], values[1], values[2], values[3]); 
    vec4 newval = valvec * transform;
    for (int i = 0; i < 4; i++)
        values[i] = newval[i];
}

void Scene::rightmultiply(const mat4 &M, stack<mat4> &transfstack)
{
    mat4 &T = transfstack.top();
    T = M * T;
}

void Scene::readfile(string fname)
{
    string str, cmd;
    ifstream in;
    in.open(fname.c_str()); 
    if (!in.is_open())
    {
        cerr << "Unable to Open Input Data File " << fname << "\n";
        throw 2;
    }
    // I need to implement a matrix stack to store transforms.  
    // This is done using standard STL Templates 
    stack <mat4> transfstack; 
    transfstack.push(mat4(1.0));  // identity
    vec3 ambient(.2,.2,.2);
    vec3 emission(0,0,0);
    vec3 diffuse(0,0,0);
    vec3 specular(0,0,0);
    float shininess = 0;
    vec3 attenuation(1,0,0);

    getline (in, str); 
    for (; in; getline (in, str))
    {
        if ((str.find_first_not_of(" \t\r\n") == string::npos) || (str[0] == '#'))
            continue;

        stringstream s(str);
        s >> cmd;
        int i;
        float values[10]; // Position and color for light, colors for others
                          // Up to 10 params for cameras.  
        bool validinput;  // Validity of input 

        // Process the light, add it to database.
        // Lighting Command
        if (cmd == "directional" || cmd == "point")
        {
            validinput = readvals(s, 6, values); // Position/color for lts.
            if (validinput)
            {
                Light l;
                if (cmd == "directional")
                    l.ltype = Light::directional;
                else
                    l.ltype = Light::point;
                for (int i = 0; i < 3; i++)
                {
                    l.pos[i] = values[i];
                    l.color[i] = values[3 + i];
                }
                l.attenuation = attenuation;
                lights.push_back(l);
            }
        }
        else if (cmd == "attenuation")
        {
            validinput = readvals(s, 3, values);
            if (validinput)
            {
                for (i = 0; i < 3; i++)
                {
                    attenuation[i] = values[i]; 
                }
            }
        }
        else if (cmd == "ambient")
        {
            validinput = readvals(s, 3, values);
            if (validinput)
            {
                for (i = 0; i < 3; i++)
                {
                    ambient[i] = values[i]; 
                }
            }
        }
        // Material Commands
        else if (cmd == "diffuse")
        {
            validinput = readvals(s, 3, values); 
            if (validinput)
            {
                for (i = 0; i < 3; i++)
                {
                    diffuse[i] = values[i]; 
                }
            }
        }
        else if (cmd == "specular")
        {
            validinput = readvals(s, 3, values); 
            if (validinput)
            {
                for (i = 0; i < 3; i++)
                {
                    specular[i] = values[i]; 
                }
            }
        }
        else if (cmd == "emission")
        {
            validinput = readvals(s, 3, values); 
            if (validinput)
            {
                for (i = 0; i < 3; i++)
                {
                    emission[i] = values[i]; 
                }
            }
        }
        else if (cmd == "shininess")
        {
            validinput = readvals(s, 1, values); 
            if (validinput)
            {
                shininess = values[0];
            }
        }
        else if (cmd == "size")
        {
            validinput = readvals(s, 2, values); 
            if (validinput)
            {
                w = (int)values[0];
                h = (int)values[1];
            }
        }
        else if (cmd == "maxdepth")
        {
            validinput = readvals(s, 1, values); 
            if (validinput)
            {
                maxdepth = values[0];
            }
        }
        else if (cmd == "output")
        {
            string outfile;
            s >> outfile;
            if (s.fail())
            {
                cerr << "Failed reading output file, will skip\n"; 
                continue;
            }
            saveFile = outfile;
        }
        else if (cmd == "camera")
        {
            validinput = readvals(s, 10, values); // 10 values eye cen up fov
            if (validinput)
            {
                cam.eye[0] = values[0];
                cam.eye[1] = values[1];
                cam.eye[2] = values[2];
                cam.center[0] = values[3];
                cam.center[1] = values[4];
                cam.center[2] = values[5];
                cam.up = Transform::upvector(vec3(values[6],values[7],values[8]), cam.center-cam.eye);
                cam.fovy = values[9];
            }
        }
        else if (cmd == "maxverts")
        {
            if (maxverts)
            {
                cerr << "Already defined maxverts parameter. Will ignore new value\n";
                continue;
            }
            validinput = readvals(s, 1, values);
            if (validinput)
            {
                maxverts = (int)values[0];
                vertex = new vec3[maxverts];
            }
        }
        else if (cmd == "vertex")
        {
            if (!vertex)
            {
                cerr << "Vertex array was not allocated\n";
                continue;
            }
            if (usedverts == maxverts)
            {
                cerr << "Vertex array already full\n";
                continue;
            }
            validinput = readvals(s, 3, values);
            if (validinput)
            {
                vertex[usedverts].x = values[0];
                vertex[usedverts].y = values[1];
                vertex[usedverts].z = values[2];
                ++usedverts;
            }
        }
        else if (cmd == "maxvertnorms")
        {
            if (maxvertnorms)
            {
                cerr << "Already defined maxvertnorms parameter. Will ignore new value\n";
                continue;
            }
            validinput = readvals(s, 1, values);
            if (validinput)
            {
                maxvertnorms = (int)values[0];
                vertnorm = new VertNorm[maxvertnorms];
            }
        }
        else if (cmd == "vertexnormal")
        {
            if (!vertnorm)
            {
                cerr << "Vertnormal array was not allocated\n";
                continue;
            }
            if (usedvertnorms == maxvertnorms)
            {
                cerr << "Vertnorm array already full\n";
                continue;
            }
            validinput = readvals(s, 6, values);
            if (validinput)
            {
                vertnorm[usedvertnorms].v.x = values[0];
                vertnorm[usedvertnorms].v.y = values[1];
                vertnorm[usedvertnorms].v.z = values[2];
                vertnorm[usedvertnorms].n.x = values[3];
                vertnorm[usedvertnorms].n.y = values[4];
                vertnorm[usedvertnorms].n.z = values[5];
                ++usedvertnorms;
            }
        }
        else if (cmd == "tri" || cmd == "trinormal")
        {
            validinput = readvals(s, 3, values); 
            if (validinput)
            {
                Object *obj = new Object;
                if (cmd == "tri")
                {
                    obj->type = tri;
                    for (int i=0; i<3; i++)
                    {
                        if (values[i] >= usedverts)
                        {
                            cerr << "Incorrect vertex index " << values[i] << endl;
                            delete obj;
                            throw 1;
                        }
                        obj->vert[i] = vertex[(int)values[i]];
                    }
                    // Transform vertex
                    for (int i=0; i<3; i++)
                    {
                        vec4 tmp(obj->vert[i], 1.0);
                        tmp = tmp * transfstack.top();
                        obj->vert[i][0] = tmp[0]/tmp[3];
                        obj->vert[i][1] = tmp[1]/tmp[3];
                        obj->vert[i][2] = tmp[2]/tmp[3];
                    }
                    // Calculate normal
                    obj->vert[3] = glm::normalize(glm::cross(obj->vert[1]-obj->vert[0], obj->vert[2]-obj->vert[0]));
                }
                else if (cmd == "trinormal")
                {
                    obj->type = trinormal;
                    for (int i=0; i<3; i++)
                    {
                        if (values[i] >= usedvertnorms)
                        {
                            cerr << "Incorrect vertexnorm index " << values[i] << endl;
                            delete obj;
                            throw 1;
                        }
                        obj->vert_norm[i] = vertnorm[(int)values[i]];
                    }
                    // Transform vertnorm
                    for (int i=0; i<3; i++)
                    {
                        vec4 tmp(obj->vert_norm[i].v, 1.0);
                        tmp = tmp * transfstack.top();
                        obj->vert_norm[i].v[0] = tmp[0]/tmp[3];
                        obj->vert_norm[i].v[1] = tmp[1]/tmp[3];
                        obj->vert_norm[i].v[2] = tmp[2]/tmp[3];
                        tmp = vec4(obj->vert_norm[i].n, 0.0);
                        tmp = glm::normalize(tmp * glm::inverseTranspose(transfstack.top()));
                        obj->vert_norm[i].n[0] = tmp[0];
                        obj->vert_norm[i].n[1] = tmp[1];
                        obj->vert_norm[i].n[2] = tmp[2];
                    }
                }
                obj->ambient = ambient;
                obj->emission = emission;
                obj->diffuse = diffuse;
                obj->specular = specular;
                obj->shininess = shininess;
                objs.push_back(obj);
            }
        }
        else if (cmd == "sphere")
        {
            validinput = readvals(s, 4, values); 
            if (validinput)
            {
                Object *obj = new Object;

                obj->type = sphere;
                for (int i = 0; i < 4; i++)
                {
                    obj->pos_r[i] = values[i];
                }
                
                // Set the object's transform
                obj->transform = transfstack.top();

                obj->ambient = ambient;
                obj->emission = emission;
                obj->diffuse = diffuse;
                obj->specular = specular;
                obj->shininess = shininess;
                objs.push_back(obj);
            }
        }
        else if (cmd == "translate")
        {
            validinput = readvals(s, 3, values); 
            if (validinput)
            {
                mat4 M = Transform::translate(values[0], values[1], values[2]);
                rightmultiply(M, transfstack);
            }
        }
        else if (cmd == "scale")
        {
            validinput = readvals(s, 3, values); 
            if (validinput)
            {
                mat4 M = Transform::scale(values[0], values[1], values[2]);
                rightmultiply(M, transfstack);
            }
        }
        else if (cmd == "rotate")
        {
            validinput = readvals(s,4,values); 
            if (validinput)
            {
                mat3 R = Transform::rotate(values[3], vec3(values[0],values[1],values[2]));
                mat4 M(1.0);
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        M[i][j] = R[i][j];
                rightmultiply(M, transfstack);
            }
        }
        else if (cmd == "pushTransform")
        {
            transfstack.push(transfstack.top()); 
        }
        else if (cmd == "popTransform")
        {
            if (transfstack.size() <= 1)
            {
                cerr << "Stack has no elements.  Cannot Pop\n"; 
            }
            else
            {
                transfstack.pop(); 
            }
        }
        else
        {
            cerr << "Unknown Command: " << cmd << " Skipping \n"; 
        }
    }
}

bool Scene::readvals(stringstream &s, const int numvals, float *values)
{
    for (int i = 0; i < numvals; i++)
    {
        s >> values[i]; 
        if (s.fail())
        {
            cout << "Failed reading value " << i << " will skip\n"; 
            return false;
        }
    }
    return true; 
}

void Scene::createBV(void )
{
    // create BV for every object and store in bvs
    for (int i = 0, n = objs.size(); i < n; i++)
    {
        Object &obj = *objs[i];
        obj.bv = BVp(new BoundingVolume);
        obj.bv->obj = &obj;
        switch (obj.type)
        {
            case tri:
            {
                vec3 min = obj.vert[0];
                vec3 max = obj.vert[0];
                for (int j = 1; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                    {
                        if (min[k] > obj.vert[j][k])
                            min[k] = obj.vert[j][k];
                        if (max[k] < obj.vert[j][k])
                            max[k] = obj.vert[j][k];
                    }
                obj.bv->lower = min;
                obj.bv->upper = max;
                break;
            }
            case trinormal:
            {
                vec3 min = obj.vert_norm[0].v;
                vec3 max = obj.vert_norm[0].v;
                for (int j = 1; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                    {
                        if (min[k] > obj.vert_norm[j].v[k])
                            min[k] = obj.vert_norm[j].v[k];
                        if (max[k] < obj.vert_norm[j].v[k])
                            max[k] = obj.vert_norm[j].v[k];
                    }
                obj.bv->lower = min;
                obj.bv->upper = max;
                break;
            }
            case sphere:
            {
                vec4 corn[8];
                corn[0][0] = obj.pos_r[0] - obj.pos_r[3];
                corn[0][1] = obj.pos_r[1] - obj.pos_r[3];
                corn[0][2] = obj.pos_r[2] - obj.pos_r[3];
                corn[0][3] = 1.f;
                
                corn[7][0] = obj.pos_r[0] + obj.pos_r[3];
                corn[7][1] = obj.pos_r[1] + obj.pos_r[3];
                corn[7][2] = obj.pos_r[2] + obj.pos_r[3];
                corn[7][3] = 1.f;
                
                corn[1] = corn[0];
                corn[1][0] += 2.f*obj.pos_r[3];
                
                corn[2] = corn[0];
                corn[2][2] += 2.f*obj.pos_r[3];
                
                corn[4] = corn[0];
                corn[4][1] += 2.f*obj.pos_r[3];
                
                corn[3] = corn[7];
                corn[3][1] -= 2.f*obj.pos_r[3];
                
                corn[5] = corn[7];
                corn[5][2] -= 2.f*obj.pos_r[3];
                
                corn[6] = corn[7];
                corn[6][0] -= 2.f*obj.pos_r[3];
                
                for (int j = 0; j < 8; j++)
                {
                    corn[j] = corn[j] * obj.transform;
                    corn[j][0] /= corn[j][3];
                    corn[j][1] /= corn[j][3];
                    corn[j][2] /= corn[j][3];
                }
                
                obj.bv->lower = vec3(corn[0][0], corn[0][1], corn[0][2]);
                obj.bv->upper = vec3(corn[0][0], corn[0][1], corn[0][2]);
                for (int j = 1; j < 8; j++)
                    for (int k = 0; k < 3; k++)
                    {
                        if (obj.bv->lower[k] > corn[j][k])
                            obj.bv->lower[k] = corn[j][k];
                        if (obj.bv->upper[k] < corn[j][k])
                            obj.bv->upper[k] = corn[j][k];
                    }
                break;
            }
        }
        bvs.push_back(obj.bv);
    }
    
    // build BV tree
    rootBV = buildBvTree(0, bvs.size(), 0);
    
    // DEBUG print tree
//     vector<BVp> stack;
//     stack.push_back(rootBV);
//     while (!stack.empty())
//     {
//         BVp cur = stack.back();
//         stack.pop_back();
//         if (cur->child2)
//             stack.push_back(cur->child2);
//         if (cur->child1)
//             stack.push_back(cur->child1);
//         printf("\nnode %p: %p %p\n", cur.get(), cur->child1.get(), cur->child2.get());
//     }
}

BVp Scene::buildBvTree(int start, int end, int axis)
{
    deque< BVp >::iterator s = bvs.begin();
    advance(s, start);
    if (end - start == 1)
        return *s;
    deque< BVp >::iterator e = bvs.begin();
    advance(e, end);
    BVp node = BVp(new BoundingVolume);
    for (deque< BVp >::iterator it = s; it != e; ++it)
        node->expand(*it);
    BoundingVolume::sort(bvs, start, end, axis);
    axis = (axis+1)%3;
    node->addChild(buildBvTree(start, (start + end)/2, axis));
    node->addChild(buildBvTree((start + end)/2, end, axis));
    return node;
}
