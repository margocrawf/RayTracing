#include "vec2.h"
#include "vec3.h"
#include "vec4.h"
#include "mat4x4.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include <algorithm>


class Material
{
    vec3 color;
public:
    Material(vec3 color) : color(color) {
    }

    virtual vec3 shade(vec3 position, vec3 normal, vec3 viewDir,
                       vec3 lightDir, vec3 powerDensity) {
        return getColor(position, normal, viewDir);
    }

	virtual vec3 getColor(
		vec3 position,
		vec3 normal,
		vec3 viewDir)
	{
		return color * std::max(normal.dot(viewDir), 0.0f); //multiply by dot product
	}
};

class DiffuseMaterial: public Material
{
    vec3 kd; // diffuse light constant

public:
    DiffuseMaterial(vec3 color) : 
        Material(color) {
        kd = vec3(1,1,0);
    }

    vec3 shade(vec3 position, vec3 normal, vec3 viewDir,
               vec3 lightDir, vec3 powerDensity) {
        // L = powerDensity . k_d(n * lightDir)+
        return powerDensity * ( kd * -(normal.dot(lightDir)));
    }
};

/*
class Metal: public Material
{
};
*/

class Wood : public Material
{
	float scale;
	float turbulence;
	float period;
	float sharpness;
public:
	Wood():
		Material(vec3(1, 1, 1))
	{
		scale = 16;
		turbulence = 500;
		period = 8;
		sharpness = 10;
	}

    float snoise(vec3 r) {
        unsigned int x = 0x0625DF73;
        unsigned int y = 0xD1B84B45;
        unsigned int z = 0x152AD8D0;
        float f = 0;
        for(int i=0; i<32; i++) {
            vec3 s(	x/(float)0xffffffff,
                    y/(float)0xffffffff, 
                    z/(float)0xffffffff);
            f += sin(s.dot(r));
            x = x << 1 | x >> 31;
            y = y << 1 | y >> 31;
            z = z << 1 | z >> 31;
        }
        return f / 64.0 + 0.5;
    }

	virtual vec3 getColor(
		vec3 position,
		vec3 normal,
		vec3 viewDir)
	{
		//return normal;
		float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence + 10000.0;
		w -= int(w);
		return (vec3(1, 0.3, 0) * w + vec3(0.35, 0.1, 0.05) * (1-w)) * normal.dot(viewDir);
	}
};
