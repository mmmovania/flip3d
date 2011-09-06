/*
 *  corrector.h
 *  flip3D
 */

#include "common.h"
#include "sorter.h"
#include <vector>

namespace corrector {
	void resample( sorter *sort, FLOAT p[3], FLOAT u[3], FLOAT re );
	void correct( sorter *sort, std::vector<particle *> &particle, FLOAT dt, FLOAT re);
};