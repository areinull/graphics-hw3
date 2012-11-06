#ifndef _SCENE_H_
#define _SCENE_H_

#include "Transform.h"
#include "Camera.h"
#include "Ray.h"
#include "Intersection.h"
#include "Light.h"
#include <string>
#include <stack>
#include <vector>
#include <sstream>
#include <FreeImage.h>
using namespace std;

class Scene
{
    enum Shape {tri, trinormal, sphere};
    struct VertNorm
    {
        vec3 v;
        vec3 n;
    };
    struct Object
    {
        Shape type;
        vec4 pos_r; // position and radius for sphere
        vec3 vert[4]; // 3 vertex + interpolated normal for triangle
        VertNorm vert_norm[3]; // 3 vertex + 3 normals for trinormal
        vec3 ambient;
        vec3 emission;
        vec3 diffuse;
        vec3 specular;
        float shininess;
        mat4 transform;
    };
    
    // Screen size
    int w, h;
    float aspect;
    int maxdepth;
    BYTE *fb; // framebuffer
    Camera cam; // camera
    string saveFile;
    int maxverts;
    int usedverts;
    vec3 *vertex; // vertex array
    int maxvertnorms;
    int usedvertnorms;
    VertNorm *vertnorm;
    vector<Object*> objs; // objects array
    vector<Light> lights; // lights array
    
    Scene(void);
    Scene(const Scene&);
    
    void readfile(string fname);
    void matransform(stack<mat4> &transfstack, float *values);
    void rightmultiply(const mat4 &M, stack<mat4> &transfstack);
    bool readvals(stringstream &s, const int numvals, float *values);
    Ray getRay(int i, int j);
    Intersection intersect(Ray& r);
    vec3 findColor(const Intersection& hit, int depth);

public:
    Scene(string fname);
    ~Scene(void);
    
    void saveScreenshot(void);
    void render(void);
    void worker(int x);
};

#endif // _SCENE_H_
