/*
 *  utility.h
 *  flip3D
 *
 */

#include "common.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifndef _UTILITY_H
#define _UTILITY_H

#define max(i,j) (i>j?i:j)
#define min(i,j) (i>j?j:i)

#define FOR_EVERY_X_FLOW	for( int i=0; i<N+1; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N; k++ ) {
#define FOR_EVERY_Y_FLOW	for( int i=0; i<N; i++ ) for( int j=0; j<N+1; j++ ) for( int k=0; k<N; k++ ) {
#define FOR_EVERY_Z_FLOW	for( int i=0; i<N; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N+1; k++ ) {
#define FOR_EVERY_CELL(n)		for( int i=0; i<n; i++ ) for( int j=0; j<n; j++ ) for( int k=0; k<n; k++ ) {
#define END_FOR }

#ifdef _OPENMP
#include <omp.h>
#define OPENMP_FOR		_Pragma("omp parallel for" )
#define OPENMP_SECTION  _Pragma("omp section" )
#define OPENMP_BEGIN	_Pragma("omp parallel" ) {
#define OPENMP_END		}
#define OPENMP_FOR_P	_Pragma("omp for" )
#else
#define OPENMP_FOR
#define OPENMP_SECTION
#define OPENMP_BEGIN
#define OPENMP_END
#define OPENMP_FOR_P
#endif

template <class T> T *** alloc3D( int w, int h, int d ) {
	T *** field = new T **[w+1];
	for( int i=0; i<w; i++ ) {
		field[i] = new T*[h+1];
		for( int j=0; j<h; j++ ) {
			field[i][j] = new T[d];
		}
		field[i][h] = NULL;
	}
	field[w] = NULL;	
	return field;
}

template <class T> void free3D( T ***ptr ) {
	for( int i=0; ptr[i]!=NULL; i++ ) {
		for( int j=0; ptr[i][j]!=NULL; j++ ) delete [] ptr[i][j];
		delete [] ptr[i];
	}
	delete [] ptr;
}

unsigned long getMicroseconds();
double dumptime();

static void dump(const char *format, ...) {
	va_list args;
    
	FILE *console = fopen( "render/log/console.out", "a" );
	va_start(args, format);
	vfprintf(console, format, args);
	va_end(args);
	fclose(console);
    
	va_start(args, format);
	vprintf(format, args);
	va_end(args);
}

void my_rand_shuffle( std::vector<ipos> &waters );

FLOAT length( FLOAT p0[3], FLOAT p1[3] );
FLOAT length2( FLOAT p0[3], FLOAT p1[3] );
FLOAT hypot2( FLOAT a, FLOAT b, FLOAT c );

#endif