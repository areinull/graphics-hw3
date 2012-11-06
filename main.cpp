#include <iostream>
#include <FreeImage.h>
#include "Scene.h"

using namespace std;

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        cerr << "Usage: raytracer scenefile\n"; 
        exit(-1); 
    }

    FreeImage_Initialise();

    Scene sc(argv[1]);
    sc.render();
    sc.saveScreenshot();

    FreeImage_DeInitialise();
    return 0;
}
