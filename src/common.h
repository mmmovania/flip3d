/*
 *  common.h
 *  flip3D
 *
 */

#include <vector>

#ifndef _COMMON_H
#define _COMMON_H

#define FLOAT	float
#define AIR		0
#define FLUID	1
#define WALL	2

#define BOX		0
#define SPHERE	1

#define GLASS	1
#define GRAY	2
#define RED		3

#define PI          3.14159265

struct _GLMmodel;
typedef struct {
	char type;
	char shape;
	char material;
	bool visible;
	FLOAT r;
	FLOAT c[3];
	FLOAT p[2][3];
} Object;

typedef struct _particle {
	FLOAT p[3];
	FLOAT u[3];
	FLOAT n[3];
	char type;
	char visible;
	char remove;
	char thinparticle;
	FLOAT tmp[2][3];
	FLOAT m;
	FLOAT dens;
} particle;

typedef struct _particle_iterator {
	particle *p;
	_particle_iterator *next;
} particle_iterator;

typedef struct _ipos {
	int i; int j; int k;
} ipos;

#endif
