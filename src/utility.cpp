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
#include <sys/time.h>

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

unsigned long getMicroseconds() {
	struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

double dumptime() {
	static unsigned prevMicroSec = getMicroseconds();
	unsigned curMicroSec = getMicroseconds();
	double res = (curMicroSec - prevMicroSec)/1000000.0;
	prevMicroSec = curMicroSec;
	return res;
}