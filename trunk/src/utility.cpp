/*
 *  utility.cpp
 *  flip3D
 *
 */

#include "utility.h"
#undef min
#undef max
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FLOAT length( FLOAT p0[3], FLOAT p1[3] ) {
	return hypotf(hypotf(p0[0]-p1[0],p0[1]-p1[1]),p0[2]-p1[2]);
}

FLOAT hypot2( FLOAT a, FLOAT b, FLOAT c ) {
    return a*a + b*b + c*c;
}

FLOAT length2( FLOAT p0[3], FLOAT p1[3] ) {
    return hypot2(p0[0]-p1[0],p0[1]-p1[1],p0[2]-p1[2]);
}

void my_rand_shuffle( std::vector<ipos> &waters ) {
	random_shuffle( waters.begin(), waters.end() );
}