/*
 *  exporter.h
 *  mdflip
 *
 *  Created by mscp on 1/9/11.
 *  Copyright 2011 さf. All rights reserved.
 *
 */

#include <vector>
#include "common.h"

namespace exporter {
	void write3D( int step, std::vector<double> &vertices, std::vector<double> &normals, std::vector<int> &faces, 
				 std::vector<Object> &objects, std::vector<particle *> &particles, int gn, FLOAT density );
}