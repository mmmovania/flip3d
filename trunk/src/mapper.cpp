/*
 *  mapper.cpp
 *  flip3D
 */

#include "mapper.h"
#include "utility.h"
#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

#define RE			1.4
#define FOR_EVERY_PARTICLE for( int n=0; n<particles.size(); n++ ) { particle *p = particles[n];

void mapper::mapP2G( sorter *sort, vector<particle *> &particles, FLOAT ****grid, int gn ) {
	
	// Compute Mapping
	OPENMP_FOR FOR_EVERY_CELL(gn+1) {
		
		// Variales for Particle Sorter
		vector<particle *> neighbors;
		
		// Map X Grids
		if( j < gn && k < gn) {
			FLOAT px[3] = { i, j+0.5, k+0.5 };
			FLOAT sumw = 0.0;
			FLOAT sumx = 0.0;
			neighbors = sort->getNeigboringParticles_wall(i,j,k,1,2,2);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *p = neighbors[n];
				if( p->type == FLUID ) {
					FLOAT x = fmax(0,fmin(gn,gn*p->p[0]));
					FLOAT y = fmax(0,fmin(gn,gn*p->p[1]));
					FLOAT z = fmax(0,fmin(gn,gn*p->p[2]));
					FLOAT pos[3] = { x, y, z };
					FLOAT w = p->m * kernel::sharp_kernel(length2(pos,px),RE);
					sumx += w*p->u[0];
					sumw += w;
				}
			}
			grid[0][i][j][k] = sumw ? sumx/sumw : 0.0;
		}
		
		// Map Y Grids
		if( i < gn && k < gn ) {
			FLOAT py[3] = { i+0.5, j, k+0.5 };
			FLOAT sumw = 0.0;
			FLOAT sumy = 0.0;
			neighbors = sort->getNeigboringParticles_wall(i,j,k,2,1,2);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *p = neighbors[n];
				if( p->type == FLUID ) {
					FLOAT x = fmax(0,fmin(gn,gn*p->p[0]));
					FLOAT y = fmax(0,fmin(gn,gn*p->p[1]));
					FLOAT z = fmax(0,fmin(gn,gn*p->p[2]));
					FLOAT pos[3] = { x, y, z };
					FLOAT w = p->m * kernel::sharp_kernel(length2(pos,py),RE);
					sumy += w*p->u[1];
					sumw += w;
				}
			}
			grid[1][i][j][k] = sumw ? sumy/sumw : 0.0;
		}
		
		// Map Z Grids
		if( i < gn && j < gn ) {
			FLOAT pz[3] = { i+0.5, j+0.5, k };
			FLOAT sumw = 0.0;
			FLOAT sumz = 0.0;
			neighbors = sort->getNeigboringParticles_wall(i,j,k,2,2,1);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *p = neighbors[n];
				if( p->type == FLUID ) {
					FLOAT x = fmax(0,fmin(gn,gn*p->p[0]));
					FLOAT y = fmax(0,fmin(gn,gn*p->p[1]));
					FLOAT z = fmax(0,fmin(gn,gn*p->p[2]));
					FLOAT pos[3] = { x, y, z };
					FLOAT w = p->m * kernel::sharp_kernel(length2(pos,pz),RE);
					sumz += w*p->u[2];
					sumw += w;
				}
			}
			grid[2][i][j][k] = sumw ? sumz/sumw : 0.0;
		}
	} END_FOR
}

void mapper::mapG2P( vector<particle *> &particles, FLOAT ****grid, int gn ) {
	OPENMP_FOR FOR_EVERY_PARTICLE {
		fetchVelocity( p->p, p->u, grid, gn );
	} END_FOR;
}

static FLOAT linear ( FLOAT ***q, FLOAT x, FLOAT y, FLOAT z, int w, int h, int d ) {
	x = fmax(0.0,fmin(w,x));
	y = fmax(0.0,fmin(h,y));
	z = fmax(0.0,fmin(d,z));
	int i = min(x,w-2);
	int j = min(y,h-2);
	int k = min(z,h-2);
	
	return	(k+1-z)*(((i+1-x)*q[i][j][k]+(x-i)*q[i+1][j][k])*(j+1-y) + ((i+1-x)*q[i][j+1][k]+(x-i)*q[i+1][j+1][k])*(y-j)) +
			(z-k)*(((i+1-x)*q[i][j][k+1]+(x-i)*q[i+1][j][k+1])*(j+1-y) + ((i+1-x)*q[i][j+1][k+1]+(x-i)*q[i+1][j+1][k+1])*(y-j));
}

void mapper::fetchVelocity( FLOAT p[3], FLOAT u[3], FLOAT ****grid, int gn ) {
	u[0] = linear( grid[0], gn*p[0], gn*p[1]-0.5, gn*p[2]-0.5, gn+1, gn, gn );
	u[1] = linear( grid[1], gn*p[0]-0.5, gn*p[1], gn*p[2]-0.5, gn, gn+1, gn );
	u[2] = linear( grid[2], gn*p[0]-0.5, gn*p[1]-0.5, gn*p[2], gn, gn, gn+1 );
}

