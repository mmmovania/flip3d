/*
 *  solver.h
 *  flip3D
 *
 */

#include "common.h"

namespace solver {
	void setSubcell( char value );
	FLOAT solve( char ***A, FLOAT ***L, FLOAT ***x, FLOAT ***b, int n );
}
