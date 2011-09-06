/*
 *  corrector.cpp
 *  flip3D
 */

#include "corrector.h"
#include "kernel.h"
#include "utility.h"
#include <math.h>

using namespace std;

#define SPRING		50.0

void corrector::resample( sorter *sort, FLOAT p[3], FLOAT u[3], FLOAT re ) {
	// Variables for Neighboring Particles
	std::vector<particle *> neighbors;
	int cell_size = sort->getCellSize();
	FLOAT wsum = 0.0;
	FLOAT save[3] = { u[0], u[1], u[2] };
	u[0] = u[1] = u[2] = 0.0;
	
	// Gather Neighboring Particles
	neighbors = sort->getNeigboringParticles_cell(fmax(0,fmin(cell_size-1,cell_size*p[0])),
												  fmax(0,fmin(cell_size-1,cell_size*p[1])),
												  fmax(0,fmin(cell_size-1,cell_size*p[2])),1,1,1);
	for( int n=0; n<neighbors.size(); n++ ) {
		particle *np = neighbors[n];
		if( np->type == FLUID ) {
			FLOAT dist2 = length2(p,np->p);
			FLOAT w = np->m * kernel::sharp_kernel(dist2,re);
			u[0] += w * np->u[0];
			u[1] += w * np->u[1];
			u[2] += w * np->u[2];
			wsum += w;
		}
	}
	if( wsum ) {
		u[0] /= wsum;
		u[1] /= wsum;
		u[2] /= wsum;
	} else {
		u[0] = save[0];
		u[1] = save[1];
		u[2] = save[2];
	}
}

void corrector::correct( sorter *sort, std::vector<particle *> &particles, FLOAT dt, FLOAT re ) {
	// Variables for Neighboring Particles
	int cell_size = sort->getCellSize();
	sort->sort(particles);
	
	// Compute Pseudo Moved Point
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {		
		if( particles[n]->type == FLUID ) {
			particle *p = particles[n];
			FLOAT spring[3] = { 0.0, 0.0, 0.0 };
			FLOAT x = max(0,min(cell_size,cell_size*p->p[0]));
			FLOAT y = max(0,min(cell_size,cell_size*p->p[1]));
			FLOAT z = max(0,min(cell_size,cell_size*p->p[2]));
			std::vector<particle *> neighbors = sort->getNeigboringParticles_cell(x,y,z,1,1,1);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *np = neighbors[n];
				if( p != np ) {
					FLOAT dist = length(p->p,np->p);
					FLOAT w = SPRING * np->m * kernel::smooth_kernel(dist*dist,re);
					if( dist > 0.1*re ) {
						spring[0] += w * (p->p[0]-np->p[0]) / dist * re;
						spring[1] += w * (p->p[1]-np->p[1]) / dist * re;
						spring[2] += w * (p->p[2]-np->p[2]) / dist * re;
					} else {
						if( np->type == FLUID ) {
							spring[0] += 0.01*re/dt*(rand()%101)/100.0;
							spring[1] += 0.01*re/dt*(rand()%101)/100.0;
							spring[2] += 0.01*re/dt*(rand()%101)/100.0;
						} else {
							spring[0] += 0.05*re/dt*np->n[0];
							spring[1] += 0.05*re/dt*np->n[1];
							spring[2] += 0.05*re/dt*np->n[2];
						}
					}
				}
			}
			p->tmp[0][0] = p->p[0] + dt*spring[0];
			p->tmp[0][1] = p->p[1] + dt*spring[1];
			p->tmp[0][2] = p->p[2] + dt*spring[2];
		}
	}
	
	// Resample New Velocity
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		if( particles[n]->type == FLUID ) {
			particle *p = particles[n];
			p->tmp[1][0] = p->u[0];
			p->tmp[1][1] = p->u[1];
			p->tmp[1][2] = p->u[2];
			resample( sort, p->tmp[0], p->tmp[1], re );

		}
	}
	
	// Update
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) { 
		if( particles[n]->type == FLUID ) {
			particle *p = particles[n];
			p->p[0] = p->tmp[0][0];
			p->p[1] = p->tmp[0][1];
			p->p[2] = p->tmp[0][2];
			p->u[0] = p->tmp[1][0];
			p->u[1] = p->tmp[1][1];
			p->u[2] = p->tmp[1][2];
		}
	}
}
