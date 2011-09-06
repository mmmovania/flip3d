/*
 *  mesher.h
 *  flip3D
 */

#include "common.h"
#include "sorter.h"
#include <vector>

namespace mesher {
	void generateMesh( char ***A, sorter *sort, std::vector<particle *> &particles, FLOAT density, int mg, 
					   std::vector<double> &vertices, std::vector<double> &normals, std::vector<int> &faces );
}