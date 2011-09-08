/*
 *  flip3D.cpp
 *  flip3D
 */

#include "common.h"
#include "utility.h"
#include "flip3D.h"
#include "utility.h"
#include "solver.h"
#include "corrector.h"
#include "sorter.h"
#include "mapper.h"
#include "mesher.h"
#include "kernel.h"
#include "exporter.h"
#include "write_bmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

#define	N               32
#define WRITE_FILE      1
#define WRITE_SAVE      0
#define MAX_STEP        600

#define ALPHA           0.95
#define DT              0.6e-2
#define DENSITY         0.5
#define	GRAVITY         9.8

#define FLIP            1
#define	WALL_THICK      (1.0/N)
#define SUBCELL         1

#define DISABLE_CORRECTION		0

static int pourTime = -1;
static FLOAT pourPos[2] = { 0.0, 0.0 };
static FLOAT pourRad = 0.12;

static FLOAT *** u[3] = { NULL, NULL, NULL };
static FLOAT *** u_save[3] = { NULL, NULL, NULL };
static char  *** A = NULL;
static FLOAT *** L = NULL;
static FLOAT *** press = NULL;
static sorter *sort = NULL;
static vector<particle *> particles;
static FLOAT ***wall_normal[4] = { NULL, NULL, NULL, NULL };
static vector<double> vertices;
static vector<double> normals;
static vector<int> faces;
static vector<Object> objects;

static FLOAT max_dens = 0.0;
static int step = 0;
static int gNumStuck = 0;

static unsigned prevGlobalMicroSec = 0;
static void resetSimulationTime() {
    prevGlobalMicroSec = getMicroseconds();
}

static void writeSimulationTime() {
	unsigned curMicroSec = getMicroseconds();
	double time = (curMicroSec - prevGlobalMicroSec)/1000000.0;	
	FILE *c = fopen( "render/log/times.out", "a" );
	fprintf( c, "%d %lf\n", step, time );
	fclose(c);
}

static void compute_wall_normal() {
	// Sort Particles
	sort->sort(particles);
	
	// Compute Wall Normal
	for( int n=0; n<particles.size(); n++ ) {
		particle *p = particles[n];
		int i = fmin(N-1,fmax(0,p->p[0]*N));
		int j = fmin(N-1,fmax(0,p->p[1]*N));
		int k = fmin(N-1,fmax(0,p->p[2]*N));
		wall_normal[0][i][j][k] = wall_normal[1][i][j][k] = wall_normal[2][i][j][k] = 0.0;
		p->n[0] = p->n[1] = p->n[2] = 0.0;
		if( p->type == WALL ) {
			if( p->p[0] <= 1.1*WALL_THICK ) {
				p->n[0] = 1.0;
			} 
			if( p->p[0] >= 1.0-1.1*WALL_THICK ) {
				p->n[0] = -1.0;
			} 
			if( p->p[1] <= 1.1*WALL_THICK ) {
				p->n[1] = 1.0;
			} 
			if( p->p[1] >= 1.0-1.1*WALL_THICK ) {
				p->n[1] = -1.0;
			}
			if( p->p[2] <= 1.1*WALL_THICK ) {
				p->n[2] = 1.0;
			} 
			if( p->p[2] >= 1.0-1.1*WALL_THICK ) {
				p->n[2] = -1.0;
			}
			
			if( p->n[0] == 0.0 && p->n[1] == 0.0 && p->n[2] == 0.0 ) {
				vector<particle *> neighbors = sort->getNeigboringParticles_cell(i,j,k,3,3,3);
				for( int n=0; n<neighbors.size(); n++ ) {
					particle *np = neighbors[n];
					if( p!=np && np->type == WALL ) {
						FLOAT d = length(p->p,np->p);
						FLOAT w = 1.0/d;
						p->n[0] += w*(p->p[0]-np->p[0])/d;
						p->n[1] += w*(p->p[1]-np->p[1])/d;
						p->n[2] += w*(p->p[2]-np->p[2])/d;
					}
				}
			}
		}
		FLOAT d = hypotf(hypotf(p->n[0],p->n[1]),p->n[2]);
		if( d ) {
			p->n[0] /= d;
			p->n[1] /= d;
			p->n[2] /= d;
			wall_normal[0][i][j][k] = p->n[0];
			wall_normal[1][i][j][k] = p->n[1];
			wall_normal[2][i][j][k] = p->n[2];
		}
	}
	
	sort->sort(particles);
	sort->markWater(A,wall_normal[3],DENSITY);
	
	// Compute Perimeter Normal
	FOR_EVERY_CELL(N) {
		wall_normal[3][i][j][k] = 0.0;
		if( A[i][j][k] != WALL ) {
			// For Every Nearby Cells
			int sum = 0;
			FLOAT norm[3] = { 0.0, 0.0, 0.0 };
			int neighbors[][3] = { {i-1,j,k}, {i+1,j,k}, {i,j-1,k}, {i,j+1,k}, {i,j,k-1}, {i,j,k+1} };
			for( int m=0; m<6; m++ ) {
				int si = neighbors[m][0];
				int sj = neighbors[m][1];
				int sk = neighbors[m][2];
				if( si < 0 || si > N-1 || sj < 0 || sj > N-1 || sk < 0 || sk > N-1 ) continue;
				if( A[si][sj][sk] == WALL ) {
					sum ++;
					norm[0] += wall_normal[0][si][sj][sk];
					norm[1] += wall_normal[1][si][sj][sk];
					norm[2] += wall_normal[2][si][sj][sk];
				}
			}
			if( sum > 0 ) {
				FLOAT d = hypot(hypot(norm[0],norm[1]),norm[2]);
				wall_normal[0][i][j][k] = norm[0]/d;
				wall_normal[1][i][j][k] = norm[1]/d;
				wall_normal[2][i][j][k] = norm[2]/d;
				wall_normal[3][i][j][k] = 1.0;
			}
		}
	} END_FOR;
}

static void computeDensity() {
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
	
		// Find Neighbors
		int gn = sort->getCellSize();
		if( particles[n]->type == WALL ) {
			particles[n]->dens = 1.0;
			continue;
		}
		
		FLOAT *p = particles[n]->p;
		vector<particle *> neighbors = sort->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*p[0])),
																		 fmax(0,fmin(gn-1,gn*p[1])),
																		 fmax(0,fmin(gn-1,gn*p[2])),1,1,1);
		FLOAT wsum = 0.0;
		for( int m=0; m<neighbors.size(); m++ ) {
			particle np = *neighbors[m];
			if( np.type == WALL ) continue;
			FLOAT d2 = length2(np.p,p);
			FLOAT w = np.m*kernel::smooth_kernel(d2, 4.0*DENSITY/N);
			wsum += w;
		}
		particles[n]->dens = wsum / max_dens;
	}
}

static void profileTest() {
	Object obj;
	obj.type = FLUID;
	obj.shape = BOX;
	obj.visible = true;
	obj.p[0][0] = 0.4;	obj.p[1][0] = 0.6;
	obj.p[0][1] = 0.0;	obj.p[1][1] = 0.3;
	obj.p[0][2] = 0.4;	obj.p[1][2] = 0.6;
	objects.push_back(obj);
}

static void spherePourTest() {
	Object obj;
	
	obj.type = WALL;
	obj.shape = SPHERE;
	obj.material = RED;
	obj.c[0] = 0.5;
	obj.visible = true;
	obj.c[1] = 0.4;
	obj.c[2] = 0.47;
	obj.r = 0.15;
	objects.push_back(obj);
	
	pourTime = 200;
	pourPos[0] = 0.5;
	pourPos[1] = 0.5;
	pourRad = 0.1;
	
	obj.type = FLUID;
	obj.shape = BOX;
	obj.visible = true;
	obj.p[0][0] = WALL_THICK;	obj.p[1][0] = 1.0-WALL_THICK;
	obj.p[0][1] = WALL_THICK;	obj.p[1][1] = 0.05;
	obj.p[0][2] = WALL_THICK;	obj.p[1][2] = 1.0-WALL_THICK;
	objects.push_back(obj);
}

static void waterDropTest() {
	Object obj;
	
	obj.type = FLUID;
	obj.shape = BOX;
	obj.p[0][0] = WALL_THICK;	obj.p[1][0] = 1.0-WALL_THICK;
	obj.p[0][1] = WALL_THICK;	obj.p[1][1] = 0.2;
	obj.p[0][2] = WALL_THICK;	obj.p[1][2] = 1.0-WALL_THICK;
	objects.push_back(obj);

	obj.type = FLUID;
	obj.shape = SPHERE;
	obj.c[0] = 0.5;
	obj.c[1] = 0.6;
	obj.c[2] = 0.5;
	obj.r = 0.12;
	objects.push_back(obj);
}

static void cliffPourTest() {
	Object obj;
	
	obj.type = WALL;
	obj.shape = BOX;
	obj.material = RED;
	obj.visible = true;
	obj.p[0][0] = 0.0;	obj.p[1][0] = 0.35;
	obj.p[0][1] = 0.0;	obj.p[1][1] = 0.4;
	obj.p[0][2] = -0.002;	obj.p[1][2] = 1.002;
	
	objects.push_back(obj);
	
	obj.type = FLUID;
	obj.shape = BOX;
	obj.visible = true;
	obj.p[0][0] = 0.35;			obj.p[1][0] = 1.0-WALL_THICK;
	obj.p[0][1] = WALL_THICK;	obj.p[1][1] = 0.05;
	obj.p[0][2] = WALL_THICK;	obj.p[1][2] = 1.0-WALL_THICK;
	
	objects.push_back(obj);
	
	pourTime = 200;
	pourPos[0] = 0.2;
	pourPos[1] = 0.2;
	pourRad = 0.1;
}

static void damBreakTest() {
	Object obj;
	
	obj.type = FLUID;
	obj.shape = BOX;
	obj.visible = true;
	obj.p[0][0] = 0.2;	obj.p[1][0] = 0.4;
	obj.p[0][1] = WALL_THICK;	obj.p[1][1] = 0.4;
	obj.p[0][2] = 0.2;	obj.p[1][2] = 0.8;
	
	objects.push_back(obj);
	
	obj.type = FLUID;
	obj.shape = BOX;
	obj.visible = true;
	obj.p[0][0] = WALL_THICK;	obj.p[1][0] = 1.0-WALL_THICK;
	obj.p[0][1] = WALL_THICK;	obj.p[1][1] = 0.06;
	obj.p[0][2] = WALL_THICK;	obj.p[1][2] = 1.0-WALL_THICK;
	
	objects.push_back(obj);
}


static void placeWalls() {
	Object obj;
	
	// Left Wall
	obj.type = WALL;
	obj.shape = BOX;
	obj.material = GLASS;
	obj.visible = 0;
	obj.p[0][0] = 0.0;			obj.p[1][0] = WALL_THICK;
	obj.p[0][1] = 0.0;			obj.p[1][1] = 1.0;
	obj.p[0][2] = 0.0;			obj.p[1][2] = 1.0;
	objects.push_back(obj);
	
	// Right Wall
	obj.type = WALL;
	obj.shape = BOX;
	obj.material = GLASS;
	obj.visible = 0;
	obj.p[0][0] = 1.0-WALL_THICK;	obj.p[1][0] = 1.0;
	obj.p[0][1] = 0.0;				obj.p[1][1] = 1.0;
	obj.p[0][2] = 0.0;				obj.p[1][2] = 1.0;
	objects.push_back(obj);
	
	// Floor Wall
	obj.type = WALL;
	obj.shape = BOX;
	obj.material = GRAY;
	obj.visible = 0;
	obj.p[0][0] = 0.0;	obj.p[1][0] = 1.0;
	obj.p[0][1] = 0.0;	obj.p[1][1] = WALL_THICK;
	obj.p[0][2] = 0.0;	obj.p[1][2] = 1.0;
	objects.push_back(obj);
	
	// Ceiling Wall
	obj.type = WALL;
	obj.shape = BOX;
	obj.material = GLASS;
	obj.visible = 0;
	obj.p[0][0] = 0.0;				obj.p[1][0] = 1.0;
	obj.p[0][1] = 1.0-WALL_THICK;	obj.p[1][1] = 1.0;
	obj.p[0][2] = 0.0;				obj.p[1][2] = 1.0;
	objects.push_back(obj);
	
	// Front Wall
	obj.type = WALL;
	obj.shape = BOX;
	obj.material = GLASS;
	obj.visible = 0;
	obj.p[0][0] = 0.0;	obj.p[1][0] = 1.0;
	obj.p[0][1] = 0.0;	obj.p[1][1] = 1.0;
	obj.p[0][2] = 0.0;	obj.p[1][2] = WALL_THICK;
	objects.push_back(obj);
	
	// Back Wall
	obj.type = WALL;
	obj.shape = BOX;
	obj.material = GLASS;
	obj.visible = 0;
	obj.p[0][0] = 0.0;				obj.p[1][0] = 1.0;
	obj.p[0][1] = 0.0;				obj.p[1][1] = 1.0;
	obj.p[0][2] = 1.0-WALL_THICK;	obj.p[1][2] = 1.0;
	objects.push_back(obj);
}

static void placeObjects() {
	// Place Object Wall
	placeWalls();
	
	// profileTest();
    // waterDropTest();
	// cliffPourTest();
	damBreakTest();
	// spherePourTest();
}

static void pourWater( int limit ) {
    if( step > limit ) return;
    
    int cnt = 0;
	double w = DENSITY/N;
    for( FLOAT x=w+w/2.0; x < 1.0-w/2.0; x += w ) {
         for( FLOAT z=w+w/2.0; z < 1.0-w/2.0; z += w ) {
             if( hypot(x-pourPos[0],z-pourPos[1]) < pourRad ) {
                 particle *p = new particle;
                 p->p[0] = x;
                 p->p[1] = 1.0 - WALL_THICK - 2.5*DENSITY/N;
                 p->p[2] = z;
                 p->u[0] = 0.0;
                 p->u[1] = -0.5*DENSITY/N/DT;
                 p->u[2] = 0.0;
                 p->n[0] = 0.0;
                 p->n[1] = 0.0;
                 p->n[2] = 0.0;
				 p->thinparticle = 0;
                 p->type = FLUID;
                 p->dens = max_dens;
                 p->m = 1.0;
                 particles.push_back(p);
                 cnt ++;
             }
         }
    }
}


static void reposition( sorter *sort, char ***A, vector<int> indices, vector<particle *> &particles ) {
	if( indices.empty() ) return;
	int gn = sort->getCellSize();
	
	// First Search for Deep Water
	vector<ipos> waters;
	while( waters.size() < indices.size() ) {
		FOR_EVERY_CELL(gn) {
#if 1
			if( i > 0 && A[i-1][j][k] != FLUID ) continue;
			if( i < N-1 && A[i+1][j][k] != FLUID ) continue;
			if( j > 0 && A[i][j-1][k] != FLUID ) continue;
			if( j < N-1 && A[i][j+1][k] != FLUID ) continue;
			if( k > 0 && A[i][j][k-1] != FLUID ) continue;
			if( k < N-1 && A[i][j][k+1] != FLUID ) continue;
#endif
			if( A[i][j][k] != FLUID ) continue;
			
			ipos aPos = { i, j, k };
			waters.push_back(aPos);
			if( waters.size() >= indices.size() ) {
				i = N; j = N; k = N;
			}
		} END_FOR;
		if( waters.empty() ) return;
	}
	
	// Shuffle
	my_rand_shuffle(waters);
	
	FLOAT h = 1.0/gn;
	for( int n=0; n<indices.size(); n++ ) {
		particle &p = *particles[indices[n]];
		p.p[0] = h*(waters[n].i+0.25+0.5*(rand()%101)/100);
		p.p[1] = h*(waters[n].j+0.25+0.5*(rand()%101)/100);
		p.p[2] = h*(waters[n].k+0.25+0.5*(rand()%101)/100);
	}
	
	sort->sort(particles);
	
	for( int n=0; n<indices.size(); n++ ) {
		particle &p = *particles[indices[n]];
		FLOAT u[3] = { 0.0, 0.0, 0.0 };
		corrector::resample( sort, p.p, u, h );
		for( int c=0; c<3; c++ ) p.u[c] = u[c];
	}
}

static void saveState() {
	static char path[256];
	static char command[256];
	static bool firstTime = true;
	if( firstTime ) {
		// Make Directory
		system("mkdir render/save" );
	}
	
	// Write Particles
	sprintf(path, "render/save/%d_scene_tmp.dat", step );
	
	dump( "Saving Status...");
	
	FILE *fp = fopen( path, "wb" );
    if( ! fp ) {
        dump( "Couldn't write to %s !\n", path );
        return;
    }
	int nparticle = particles.size();
	fwrite( &nparticle, sizeof(int), 1, fp );
	fwrite( &max_dens, sizeof(FLOAT), 1, fp );
	fwrite( &pourTime, sizeof(int), 1, fp );
	fwrite( pourPos, 2*sizeof(FLOAT), 1, fp );
	fwrite( &pourRad, sizeof(FLOAT), 1, fp );
	
	for( int n=0; n < particles.size(); n++ ) {
		if( n % (particles.size()/10) == 0 ) dump( ">" );
		
		// Write Particle Structure
		fwrite( particles[n], sizeof(particle), 1, fp );
	}
    
    FOR_EVERY_CELL(N) {
        // Write Pressure
        fwrite( &press[i][j][k], sizeof(FLOAT), 1, fp );
    } END_FOR;
	
	// Save Cell Normal Walls
	FOR_EVERY_CELL(N) {
		for( int n=0; n<4; n++ ) fwrite( &wall_normal[n][i][j][k], sizeof(FLOAT), 1, fp );
	} END_FOR;

	fclose(fp);
	
	if( firstTime ) {
		// Write Objects
		FILE *fp = fopen( "render/save/objects.dat", "w" );
		fprintf( fp, "%d\n", (int)objects.size() );
		for( int n=0; n < objects.size(); n++ ) {
			Object &obj = objects[n];
			fprintf( fp, "%d %d %d %d\n", obj.type, obj.shape, obj.material, obj.visible );
			fprintf( fp, "%f %f %f %f %f %f %f %f %f %f\n", obj.r, obj.c[0], obj.c[1], obj.c[2], 
					obj.p[0][0], obj.p[0][1], obj.p[0][2], obj.p[1][0], obj.p[1][1], obj.p[1][2] );
		}
		fclose(fp);
	}
	firstTime = false;
	
	sprintf(command, "mv render/save/%d_scene_tmp.dat render/save/%d_scene.dat", step, step );
	system(command);
	
	dump( " %.2f sec\n", dumptime() );
}

static void loadState( int load_step ) {
	static char path[256];
	static bool firstTime = true;
	if( firstTime ) {
		// Make Directory
		system("mkdir render/save" );
	}
	
	// Load Particles
	sprintf(path, "render/save/%d_scene.dat", load_step );
	
	FILE *fp = fopen( path, "rb" );
    if( ! fp ) {
        dump( "Couldn't read file %s, exiting...\n", path );
        exit(0);
    }
	int nparticle = 0;
	fread( &nparticle, sizeof(int), 1, fp );
	fread( &max_dens, sizeof(FLOAT), 1, fp );
	fread( &pourTime, sizeof(int), 1, fp );
	fread( pourPos, 2*sizeof(FLOAT), 1, fp );
	fread( &pourRad, sizeof(FLOAT), 1, fp );
	particles.resize(nparticle);
	for( int n=0; n < nparticle; n++ ) {
		particle *newp = new particle;
    
		// Load Basic Structure
		fread( newp, sizeof(particle), 1, fp );
        particles[n] = newp;
	}
    
    FOR_EVERY_CELL(N) {
        // Load Pressure
        fread( &press[i][j][k], sizeof(FLOAT), 1, fp );
    } END_FOR;
	
	// Save Cell Normal Walls
	FOR_EVERY_CELL(N) {
		for( int n=0; n<4; n++ ) fread( &wall_normal[n][i][j][k], sizeof(FLOAT), 1, fp );
	} END_FOR;
	
	fclose(fp);

	if( firstTime ) {
		// Write Objects
		FILE *fp = fopen( "render/save/objects.dat", "r" );
		int nobjects = 0;
		fscanf( fp, "%d\n", &nobjects );
		objects.resize(nobjects);
		for( int n=0; n < nobjects; n++ ) {
			Object obj;
			fscanf( fp, "%d %d %d %d\n", &obj.type, &obj.shape, &obj.material, &obj.visible );
			fscanf( fp, "%f %f %f %f %f %f %f %f %f %f\n", &obj.r, &obj.c[0], &obj.c[1], &obj.c[2], 
					&obj.p[0][0], &obj.p[0][1], &obj.p[0][2], &obj.p[1][0], &obj.p[1][1], &obj.p[1][2] );
			objects[n] = obj;
		}
		fclose(fp);
	}
	firstTime = false;
	step = load_step;
}

static void pushParticle( double x, double y, double z, char type ) {
	Object *inside_obj = NULL;
	for( int n=0; n<objects.size(); n++ ) {
		Object &obj = objects[n];
		
		bool found = false;
		FLOAT thickness = 3.0/N;
		if( obj.shape == BOX ) {
			if( x > obj.p[0][0] && x < obj.p[1][0] &&
			   y > obj.p[0][1] && y < obj.p[1][1] &&
			   z > obj.p[0][2] && z < obj.p[1][2] ) {
				
				if(	 obj.type == WALL &&
				   x > obj.p[0][0]+thickness && x < obj.p[1][0]-thickness &&
				   y > obj.p[0][1]+thickness && y < obj.p[1][1]-thickness &&
				   z > obj.p[0][2]+thickness && z < obj.p[1][2]-thickness ) {
					// Do nothing. Because It's too deep
					inside_obj = NULL;
					break;
				} else {
					found = true;
				}
			}
		} else if( obj.shape == SPHERE ) {
			FLOAT p[3] = { x, y, z };
			FLOAT c[3] = { obj.c[0], obj.c[1], obj.c[2] };
			FLOAT len = length(p, c);
			if( len < obj.r ) {
				if( obj.type == WALL ) {
					found = true;
					if( len < obj.r-thickness ) {
						// Do nothing. Because It's too deep
						inside_obj = NULL;
						break;
					}
				} else if( obj.type == FLUID ) {
					found = true;
				}
			}
		}
		
		if( found ) {
			if(	objects[n].type == type ) {
				inside_obj = &objects[n]; // Found
				break;
			}
		}
	}
	
	if( inside_obj ) {
		particle *p = new particle;
		p->p[0] = x + 0.01*(inside_obj->type==FLUID)*0.2*((rand()%101)/50.0-1.0)/N;
		p->p[1] = y + 0.01*(inside_obj->type==FLUID)*0.2*((rand()%101)/50.0-1.0)/N;
		p->p[2] = z + 0.01*(inside_obj->type==FLUID)*0.2*((rand()%101)/50.0-1.0)/N;
        p->u[0] = 0.0;
        p->u[1] = 0.0;
        p->u[2] = 0.0;
		p->n[0] = 0.0;
		p->n[1] = 0.0;
		p->n[2] = 0.0;
		p->thinparticle = 0;
		p->dens = 10.0;
		p->type = inside_obj->type;
		p->visible = inside_obj->visible;
		p->m = 1.0;
		particles.push_back(p);
	}
}

void flip3D::init( int load_step ) {
	// Allocate Memory
	u[0] = alloc3D<FLOAT>(N+1,N,N);
	u[1] = alloc3D<FLOAT>(N,N+1,N);
	u[2] = alloc3D<FLOAT>(N,N,N+1);
	u_save[0] = alloc3D<FLOAT>(N+1,N,N);
	u_save[1] = alloc3D<FLOAT>(N,N+1,N);
	u_save[2] = alloc3D<FLOAT>(N,N,N+1);
	A = alloc3D<char>(N,N,N);
	L = alloc3D<FLOAT>(N,N,N);
    press = alloc3D<FLOAT>(N,N,N);
	wall_normal[0] = alloc3D<FLOAT>(N,N,N);
	wall_normal[1] = alloc3D<FLOAT>(N,N,N);
	wall_normal[2] = alloc3D<FLOAT>(N,N,N);
	wall_normal[3] = alloc3D<FLOAT>(N,N,N);

	FOR_EVERY_X_FLOW {
		u_save[0][i][j][k] = u[0][i][j][k] = 0.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		u_save[1][i][j][k] = u[1][i][j][k] = 0.0;
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		u_save[2][i][j][k] = u[2][i][j][k] = 0.0;
	} END_FOR
	
	FOR_EVERY_CELL(N) {
		A[i][j][k] = AIR;
        press[i][j][k] = 0.0;
	} END_FOR
	
	// Allocate Sorter
	if( ! sort ) sort = new sorter(N);
	
	if( load_step == 0 ) {
		// Define Objects
		placeObjects();
		
		// This Is A Test Part. We Generate Pseudo Particles To Measure Maximum Particle Density
		FLOAT h = DENSITY/N;
		FOR_EVERY_CELL(10) {
			particle *p = new particle;
			p->p[0] = (i+0.5)*h;
			p->p[1] = (j+0.5)*h;
			p->p[2] = (k+0.5)*h;
			p->type = FLUID;
			p->m = 1.0;
			particles.push_back(p);
		} END_FOR
		sort->sort(particles);
		max_dens = 1.0;
		computeDensity();
		max_dens = 0.0;
		for( int n=0; n<particles.size(); n++ ) {
			particle *p = particles[n];
			max_dens = fmax(max_dens,p->dens);
			delete p;
		}
		particles.clear();
		
		// Place Fluid Particles And Walls
		double w = DENSITY*WALL_THICK;
		for( int i=0; i < N/DENSITY; i++ ) {
			for( int j=0; j < N/DENSITY; j++ ) {
				for( int k=0; k < N/DENSITY; k++ ) {
					double x = i*w+w/2.0;
					double y = j*w+w/2.0;
					double z = k*w+w/2.0;
			
					if( x > WALL_THICK && x < 1.0-WALL_THICK &&
						y > WALL_THICK && y < 1.0-WALL_THICK &&
						z > WALL_THICK && z < 1.0-WALL_THICK ) {
							pushParticle( x, y, z, FLUID );
					}					   
				}
			}
		}
		
		// Place Wall Particles And Walls
		w = 1.0/N;
		for( int i=0; i < N; i++ ) {
			for( int j=0; j < N; j++ ) {
				for( int k=0; k < N; k++ ) {
					double x = i*w+w/2.0;
					double y = j*w+w/2.0;
					double z = k*w+w/2.0;
					pushParticle( x, y, z, WALL );
				}
			}
		}
		
		// Remove Particles That Stuck On Wal Cells
		sort->sort(particles);
		sort->markWater(A,wall_normal[3],DENSITY);
		
		for( vector<particle *>::iterator iter=particles.begin(); iter!=particles.end(); ) {
			particle &p = **iter;
			if( p.type == WALL ) {
				iter++;
				continue;
			}
			int i = fmin(N-1,fmax(0,p.p[0]*N));
			int j = fmin(N-1,fmax(0,p.p[1]*N));
			int k = fmin(N-1,fmax(0,p.p[2]*N));
			if( A[i][j][k] == WALL ) {
				delete *iter;
				iter = particles.erase(iter);
			} else {
				iter ++;
			}
		}
		
		// Comput Normal for Walls
		compute_wall_normal();
		
#if _OPENMP
		printf( "OpenMP Detected.\n" );
#endif
		
        system( "mkdir render" );
		system( "rm -rf render/*" );
		system( "mkdir render/preview" );
		system( "mkdir render/log" );
	} else {
		loadState(load_step);
	}
}

static inline char wall_chk( char a ) {
	return a == WALL ? 1.0 : -1.0;
}

static void enforce_boundary() {
	// Set Boundary Velocity Zero
	FOR_EVERY_X_FLOW {
		if( i==0 || i==N ) u[0][i][j][k] = 0.0;
		if( i<N && i>0 && wall_chk(A[i][j][k])*wall_chk(A[i-1][j][k]) < 0 ) {
			u[0][i][j][k] = 0.0;
		}
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		if( j==0 || j==N ) u[1][i][j][k] = 0.0;
		if( j<N && j>0 && wall_chk(A[i][j][k])*wall_chk(A[i][j-1][k]) < 0 ) {
			u[1][i][j][k] = 0.0;
		}
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		if( k==0 || k==N ) u[2][i][j][k] = 0.0;
		if( k<N && k>0 && wall_chk(A[i][j][k])*wall_chk(A[i][j][k-1]) < 0 ) {
			u[2][i][j][k] = 0.0;
		}
	} END_FOR
}

static void project() {
	// Cell Width
	FLOAT h = 1.0/N;
	
	// Memory Allocation
	static FLOAT *** div = alloc3D<FLOAT>(N,N,N);
	
	// Compute Divergence
	FOR_EVERY_CELL(N) {
		if( A[i][j][k] == FLUID ) {
			div[i][j][k] = (u[0][i+1][j][k]-u[0][i][j][k]+
							u[1][i][j+1][k]-u[1][i][j][k]+
							u[2][i][j][k+1]-u[2][i][j][k]) / h;
		}
	} END_FOR;
	
	// Compute LevelSet
	FOR_EVERY_CELL(N) {
		L[i][j][k] = sort->levelset(i,j,k,wall_normal[3],DENSITY);
	} END_FOR;
	
	// Solve Pressure
	solver::setSubcell(SUBCELL);
	solver::solve( A, L, press, div, N );
	
	// Subtract Pressure Gradient
	FOR_EVERY_X_FLOW {
		if( i>0 && i<N ) {
			FLOAT pf = press[i][j][k];
			FLOAT pb = press[i-1][j][k];
			if( SUBCELL && L[i][j][k] * L[i-1][j][k] < 0.0 ) {
				pf = L[i][j][k] < 0.0 ? press[i][j][k] : L[i][j][k]/fmin(1.0e-3,L[i-1][j][k])*press[i-1][j][k];
				pb = L[i-1][j][k] < 0.0 ? press[i-1][j][k] : L[i-1][j][k]/fmin(1.0e-6,L[i][j][k])*press[i][j][k];
			}
			u[0][i][j][k] -= (pf-pb)/h;
		}
	} END_FOR;
	
	FOR_EVERY_Y_FLOW {
		if( j>0 && j<N ) {
			FLOAT pf = press[i][j][k];
			FLOAT pb = press[i][j-1][k];
			if( SUBCELL && L[i][j][k] * L[i][j-1][k] < 0.0 ) {
				pf = L[i][j][k] < 0.0 ? press[i][j][k] : L[i][j][k]/fmin(1.0e-3,L[i][j-1][k])*press[i][j-1][k];
				pb = L[i][j-1][k] < 0.0 ? press[i][j-1][k] : L[i][j-1][k]/fmin(1.0e-6,L[i][j][k])*press[i][j][k];
			}
			u[1][i][j][k] -= (pf-pb)/h;
		}
	} END_FOR;
	
	FOR_EVERY_Z_FLOW {
		if( k>0 && k<N ) {
			FLOAT pf = press[i][j][k];
			FLOAT pb = press[i][j][k-1];
			if( SUBCELL && L[i][j][k] * L[i][j][k-1] < 0.0 ) {
				pf = L[i][j][k] < 0.0 ? press[i][j][k] : L[i][j][k]/fmin(1.0e-3,L[i][j][k-1])*press[i][j][k-1];
				pb = L[i][j][k-1] < 0.0 ? press[i][j][k-1] : L[i][j][k-1]/fmin(1.0e-6,L[i][j][k])*press[i][j][k];
			}
			u[2][i][j][k] -= (pf-pb)/h;
		}
	} END_FOR;
}

static void add_ExtForce() {
	for( int n=0; n<particles.size(); n++ ) {
		// Add Gravity
		particles[n]->u[1] += -DT*GRAVITY;
	}
}

static void advect_particle() {
	// Advect Particle Through Grid
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		if( particles[n]->type == FLUID ) {
			FLOAT vel[3];
			mapper::fetchVelocity( particles[n]->p, vel, u, N );
			for( int k=0; k<3; k++ ) {
				particles[n]->p[k] += DT*vel[k];
			}
		}
	}
	
	// Sort
	sort->sort(particles);

	// Constraint Outer Wall
	for( int n=0; n<particles.size(); n++ ) {
		FLOAT r = WALL_THICK;
		for( int k=0; k<3; k++ ) {
			if( particles[n]->type == FLUID ) {
				particles[n]->p[k] = fmax(r,fmin(1.0-r,particles[n]->p[k]));
			}
		}
#if 1
		particle *p = particles[n];
		if( p->type == FLUID ) {
			int i = fmin(N-1,fmax(0,p->p[0]*N));
			int j = fmin(N-1,fmax(0,p->p[1]*N));
			int k = fmin(N-1,fmax(0,p->p[2]*N));
			vector<particle *> neighbors = sort->getNeigboringParticles_cell(i,j,k,1,1,1);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *np = neighbors[n];
				double re = 1.5*DENSITY/N;
				if( np->type == WALL ) {
					FLOAT dist = length(p->p,np->p);
					if( dist < re ) {
						FLOAT normal[3] = { np->n[0], np->n[1], np->n[2] };
						if( normal[0] == 0.0 && normal[1] == 0.0 && normal[2] == 0.0 && dist ) {
							for( int c=0; c<3; c++ ) normal[c] = (p->p[c]-np->p[c])/dist;
						}
						   
						p->p[0] += (re-dist)*normal[0];
						p->p[1] += (re-dist)*normal[1];
						p->p[2] += (re-dist)*normal[2];
						FLOAT dot = p->u[0] * normal[0] + p->u[1] * normal[1] + p->u[2] * normal[2];
						p->u[0] -= dot*normal[0];
						p->u[1] -= dot*normal[1];
						p->u[2] -= dot*normal[2];
					}
				}
			}
		}
#endif
	}
   
#if 1
    // Remove Particles That Stuck On The Up-Down Wall Cells...
    for( int n=0; n<particles.size(); n++ ) {
        particle &p = *particles[n];
		p.remove = 0;
		
		// Focus on Only Fluid Particle
		if( p.type != FLUID ) {
			continue;
		}
		
		// If Stuck On Wall Cells Just Repositoin
		if( A[(int)fmin(N-1,fmax(0,p.p[0]*N))][(int)fmin(N-1,fmax(0,p.p[1]*N))][(int)fmin(N-1,fmax(0,p.p[2]*N))] == WALL ) {
			p.remove = 1;
		}
		
        int i = fmin(N-3,fmax(2,p.p[0]*N));
        int j = fmin(N-3,fmax(2,p.p[1]*N));
        int k = fmin(N-3,fmax(2,p.p[2]*N));
        if( p.dens < 0.04 && (A[i][max(0,j-1)][k] == WALL || A[i][min(N-1,j+1)][k] == WALL) ) {
			// Put Into Reposition List
			p.remove = 1;
        }
    }
	
	// Reposition If Neccessary
	vector<int> reposition_indices;
	for( int n=0; n<particles.size(); n++ ) {
		if( particles[n]->remove ) {
			particles[n]->remove = 0;
			reposition_indices.push_back(n);
		}
	}
	// Store Stuck Particle Number
	gNumStuck = reposition_indices.size();
	reposition( sort, A, reposition_indices, particles );
#endif
}

static void save_grid() {
	FOR_EVERY_X_FLOW {
		u_save[0][i][j][k] = u[0][i][j][k];
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		u_save[1][i][j][k] = u[1][i][j][k];
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		u_save[2][i][j][k] = u[2][i][j][k];
	} END_FOR
}

static void subtract_grid() {
	FOR_EVERY_X_FLOW {
		u_save[0][i][j][k] = u[0][i][j][k] - u_save[0][i][j][k];
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		u_save[1][i][j][k] = u[1][i][j][k] - u_save[1][i][j][k];
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		u_save[2][i][j][k] = u[2][i][j][k] - u_save[2][i][j][k];
	} END_FOR
}

static void write_mesh() {
    // Generate Mesh
    dump( "[ MESH EXPORT ] Generating Meshes...");
    mesher::generateMesh( A, sort, particles, DENSITY, N+1, vertices, normals, faces );
    dump( "%.2f sec\n", dumptime() );

    // Write Into File
    dump( "[ MESH EXPORT ] Writing Meshes...");
    exporter::write3D( step, vertices, normals, faces, objects, particles, WALL_THICK );
    dump( "%.2f sec\n", dumptime() );
}

static void extrapolate_velocity() {
	// Mark Fluid Cell Face
	static char ***mark[3] = { alloc3D<char>(N+1,N,N), alloc3D<char>(N,N+1,N), alloc3D<char>(N,N,N+1) };
	static char ***wall_mark[3] = { alloc3D<char>(N+1,N,N), alloc3D<char>(N,N+1,N), alloc3D<char>(N,N,N+1) };
	
	OPENMP_FOR FOR_EVERY_X_FLOW {
		mark[0][i][j][k] = (i>0 && A[i-1][j][k]==FLUID) || (i<N && A[i][j][k]==FLUID);
		wall_mark[0][i][j][k] = (i<=0 || A[i-1][j][k]==WALL) && (i>=N || A[i][j][k]==WALL);
	} END_FOR;
	
	OPENMP_FOR FOR_EVERY_Y_FLOW {
		mark[1][i][j][k] = (j>0 && A[i][j-1][k]==FLUID) || (j<N && A[i][j][k]==FLUID);
		wall_mark[1][i][j][k] = (j<=0 || A[i][j-1][k]==WALL) && (j>=N || A[i][j][k]==WALL);
	} END_FOR;
	
	OPENMP_FOR FOR_EVERY_Z_FLOW {
		mark[2][i][j][k] = (k>0 && A[i][j][k-1]==FLUID) || (k<N && A[i][j][k]==FLUID);
		wall_mark[2][i][j][k] = (k<=0 || A[i][j][k-1]==WALL) && (k>=N || A[i][j][k]==WALL);
	} END_FOR;
	
	// Now Extrapolate
	OPENMP_FOR FOR_EVERY_CELL(N+1) {
		for( int n=0; n<3; n++ ) {
			if( n!=0 && i>N-1 ) continue;
			if( n!=1 && j>N-1 ) continue;
			if( n!=2 && k>N-1 ) continue;
			
			if( ! mark[n][i][j][k] && wall_mark[n][i][j][k] ) {
				int wsum = 0;
				FLOAT sum = 0.0;
				int q[][3] = { {i-1,j,k}, {i+1,j,k}, {i,j-1,k}, {i,j+1,k}, {i,j,k-1}, {i,j,k+1} };
				for( int qk=0; qk<6; qk++ ) {
					if( q[qk][0]>=0 && q[qk][0]<N+(n==0) && q[qk][1]>=0 && q[qk][1]<N+(n==1) && q[qk][2]>=0 && q[qk][2]<N+(n==2) ) {
						if( mark[n][q[qk][0]][q[qk][1]][q[qk][2]] ) {
							wsum ++;
							sum += u[n][q[qk][0]][q[qk][1]][q[qk][2]];
						}
					}
				}
				if( wsum ) u[n][i][j][k] = sum/wsum;
			}
		}
	} END_FOR;
}


static void solve_picflip() {
    // Map Particles Onto Grid
	sort->sort(particles);
	mapper::mapP2G(sort,particles,u,N);
	sort->markWater(A,wall_normal[3],DENSITY);
    
	// Solve Fluid Velocity On Grid
	save_grid();
	enforce_boundary();
	project();
	enforce_boundary();
	extrapolate_velocity();
	subtract_grid();
	
#if FLIP
	// Copy Current Velocity
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		for( int k=0; k<3; k++ ) {
			particles[n]->tmp[0][k] = particles[n]->u[k];
		}
	}
	
	// Map Changes Back To Particles
	mapper::mapG2P(particles,u_save,N);
	
	// Set Tmp As FLIP Velocity
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		for( int k=0; k<3; k++ ) {
			particles[n]->tmp[0][k] = particles[n]->u[k] + particles[n]->tmp[0][k];
		}
	}
	
	// Set u[] As PIC Velocity
	mapper::mapG2P(particles,u,N);
	
	// Interpolate
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		for( int k=0; k<3; k++ ) {
			particles[n]->u[k] = (1.0-ALPHA)*particles[n]->u[k] + ALPHA*particles[n]->tmp[0][k];
		}
	}
	
#else
	// Map Changes Back To Particles
	mapper::mapG2P(particles,u,N);
#endif // FLIP
}

static void buildRotationMatrix( FLOAT xrad, FLOAT yrad, FLOAT R[3][3] ) {
    FLOAT Rx[3][3] = { {1,0,0},{0,cos(xrad),-sin(xrad)},{0,sin(xrad),cos(xrad)} };
    FLOAT Ry[3][3] = { {cos(yrad),0,sin(yrad)}, {0,1,0}, {-sin(yrad),0,cos(yrad)} };
    FLOAT Rtmp[3][3];
    for( int i=0; i<3; i++ ) {
        for( int j=0; j<3; j++ ) {
            R[i][j] = Rtmp[i][j] = Rx[i][j];
        }
    }
    
    for( int i=0; i<3; i++ ) {
        for( int j=0; j<3; j++ ) {
            R[i][j] = 0.0;
            for( int k=0; k<3; k++ ) R[i][j] += Rtmp[i][k]*Ry[k][j];
        }
    }
}

static void transform( FLOAT p[3], FLOAT R[3][3] ) {
    FLOAT p0[3] = { p[0], p[1], p[2] };
    for( int i=0; i<3; i++ ) {
        p[i] = 0.0;
        for( int k=0; k<3; k++ ) {
            p[i] += R[i][k]*p0[k];
        }
    }
}

static void write_frame() {
    int width = 1024;
    int height = 700;
    FLOAT winrate = height/(FLOAT)width;
    
    static FLOAT *image = new FLOAT[width*height*4];
    static unsigned char *buffer = new unsigned char[width*height*4];
    for( int i=0; i<width; i++ ) {
        for( int j=0; j<height; j++ ) {
            for( int c=0; c<3; c++ ) {
                image[4*(i+j*width)+c] = 0.0;
            }
            image[4*(i+j*width)+3] = 1.0;
        }
    }
    
    FLOAT blueColor[] = { 0.3, 0.5, 0.8, exp(-DENSITY*N/23.0) };
    static FLOAT R[3][3];
    bool firstTime = true;
    if( firstTime ) {
        buildRotationMatrix( -0.2, 0.2, R );
        firstTime = false;
    }
    
    // Simple Point-based Rendering
    FLOAT eye = 2.0;
    FLOAT offset = 0.3;
    
    for( int n=0; n<particles.size(); n++ ) {
        if( particles[n]->type == FLUID ) {
            FLOAT p[3] = { particles[n]->p[0]-0.5, particles[n]->p[1]-0.5, particles[n]->p[2]-0.5 };
            transform(p,R); 
            FLOAT z = offset + 0.5 + p[2];
            FLOAT x = eye/(eye+z)*(p[0]+0.05);
            FLOAT y = eye/(eye+z)*(p[1]+0.25);
            if( x > -0.5 && x < 0.5 && y > -winrate/2.0 && y < winrate/2.0 ) {
                int i = width*(x+0.5);
                int j = width*(y+winrate/2.0);
                int w = 1;
                for( int ni=i-w; ni<i+1; ni++ ) {
                    for( int nj=j-w; nj<j+1; nj++ ) {
                        if( ni>0 && ni<width && nj>0 && nj<height ) {
                            for( int c=0; c<3; c++ ) {
                                image[4*(ni+nj*width)+c] = blueColor[3]*blueColor[c] + (1.0-blueColor[3])*image[4*(ni+nj*width)+c];
                            }
                        }
                    }
                }
            }
        }
    }
    
    for( int n=0; n<4*width*height; n++ ) {
        buffer[n] = 255.0*fmin(1.0,image[n]);
    }
    
    char name[256];
	sprintf( name, "render/preview/frame_%d.bmp", step );
	write_bmp( name, buffer, width, height, true );
}

void flip3D::simulateStep() {
    
	dump( "-------------- Step %d --------------\n", step++);
    
    // Pour Water
    pourWater(pourTime);
    
	// Display Env
#if _OPENMP
	static int procs = omp_get_num_procs();
	dump( "Number of threads: %d\n", procs );
#else
	dump( "OpenMP Disabled.\n" );
#endif
	
    // Write Image
    dump( "Writing Preview Image..." );
    dumptime();
	write_frame();
    dump( "%.2f sec\n", dumptime() );
    
	// Reset Simulation Time
    resetSimulationTime();
	
	// Compute Density
	dump( "Computing Density..." );
	dumptime();
	sort->sort(particles);
	computeDensity();
	dump( "%.2f sec\n", dumptime() );
	
    // Solve Fluid
	dump( "Solving Liquid Flow...");
    if( step == 1 && N > 64 ) {
        dump( "\n>>> NOTICE:\nBe advised that the first step of pressure solver really takes time.\nJust be patient :-)\n<<<\n");
    }
	add_ExtForce();
    solve_picflip();
    dump( "Took %.2f sec\n", dumptime() );
    
	// Advect Particle
	dump( "Advecting Particles...");
	advect_particle();
	dump( "%.2f sec\n", dumptime() );
	
	// Correct Position
#if ! DISABLE_CORRECTION
	dump( "Correcting Particle Position...");
	corrector::correct(sort,particles,DT,DENSITY/N);
	dump( "%.2f sec\n", dumptime() );
#endif
	
    // Write Simulation Time
    writeSimulationTime();

	// Save Current State
#if WRITE_SAVE
    saveState();
#endif
    
	// Write Mesh
#if WRITE_FILE
	write_mesh();
#endif
    
    // If Exceeds Max Step Exit
	if( step > MAX_STEP ) {
        dump( "Maximum Timestep Reached. Exiting...\n");
		exit(0);
	}
}