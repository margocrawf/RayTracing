RayCast: ray shader.h vec2.h vec3.h vec4.h mat4x4.h quadricmats.h materials.h RayCast.cpp
	g++ -std=c++11 -o RayCast RayCast.o shader.h vec2.h vec3.h vec4.h mat4x4.h quadricmats.h materials.h -lGL -lglut -lGLEW -O1

ray: RayCast.cpp
	g++ -c RayCast.cpp

