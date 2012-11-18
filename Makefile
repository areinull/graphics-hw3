CC = g++
CFLAGS = -g --std=c++0x
INCFLAGS = -I./glm-0.9.2.7 -I./include/
LDFLAGS = -lfreeimage -lboost_thread
RM = /bin/rm -f 
all: raytracer
raytracer: main.cpp Transform.cpp Transform.h Camera.cpp Camera.h Intersection.cpp Intersection.h Ray.cpp Ray.h Scene.cpp Scene.h Light.cpp Light.h BoundingVolume.cpp BoundingVolume.h
	$(CC) $(CFLAGS) -o raytracer main.cpp Transform.cpp Camera.cpp Intersection.cpp Ray.cpp Scene.cpp Light.cpp BoundingVolume.cpp $(INCFLAGS) $(LDFLAGS) 
clean: 
	$(RM) *.o raytracer *.png


 
