/*
 *  sorter.h
 *  flip3D
 */

#include "common.h"
#include <vector>
#ifndef _SORTER_H
#define _SORTER_H

class sorter {
public:
	sorter( int gn );
	~sorter();
	
	void sort( std::vector<particle *> &particles );
	std::vector<particle *> getNeigboringParticles_wall( int i, int j, int k, int w, int h, int d );
	std::vector<particle *> getNeigboringParticles_cell( int i, int j, int k, int w, int h, int d );
	FLOAT levelset( int i, int j, int k, FLOAT ***halfwall, FLOAT density );
	
	int	 getCellSize(){ return gn; }
	int	 getNumParticleAt( int i, int j, int k );
	void markWater( char ***A, FLOAT ***halfwall, FLOAT density );
	void deleteAllParticles();
	
protected:
	std::vector<particle *> ***cells;
	int gn;
};

#endif