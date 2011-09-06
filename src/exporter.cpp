/*
 *  exporter.cpp
 *  mdflip
 */

#include "exporter.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

void write_obj( int frame, vector<double> &vertices, vector<double> &normals, vector<int> &faces, FLOAT h ) {
	
	static bool firstTime = true;
	if( firstTime ) {
		// Make Directory
		system("mkdir render/obj" );
	}
	firstTime = false;
	
	char tmp[64];
	sprintf(tmp, "render/obj/%d_scene.obj", frame );
	FILE *obj_fp = fopen( tmp, "w" );
	
	// Vertices
	for( int i=0; i<vertices.size(); i+=3 ) {
		fprintf( obj_fp, "v %lf %lf %lf\n", fmin(1.0-h,fmax(h,vertices[i])), fmin(1.0-h,fmax(h,vertices[i+2])), fmin(1.0-h,fmax(h,vertices[i+1])) );
	}
	
	// Close And Open As Append
	fclose(obj_fp);
	obj_fp = fopen( tmp, "a" );
	
	// Normals
	for( int i=0; i<normals.size(); i+=3 ) {
		fprintf( obj_fp, "vn %lf %lf %lf\n", normals[i], normals[i+2], normals[i+1] );
	}
	
	// Close And Open As Append
	fclose(obj_fp);
	obj_fp = fopen( tmp, "a" );
	
	// Faces
	for( int i=0; i<faces.size(); i+=3 ) {
		fprintf( obj_fp, "f %d//%d %d//%d %d//%d\n", faces[i]+1, faces[i]+1, faces[i+1]+1, faces[i+1]+1, faces[i+2]+1, faces[i+2]+1 );
	}
    fflush(obj_fp);
	fclose(obj_fp);
}

void exporter::write3D( int step, vector<double> &vertices, vector<double> &normals, vector<int> &faces, 
					   vector<Object> &objects, std::vector<particle *> &particles, int gn, FLOAT density ) {
	write_obj( step, vertices, normals, faces, 1.0/gn );
}
