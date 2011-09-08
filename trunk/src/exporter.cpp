/*
 *  exporter.cpp
 *  mdflip
 */

#include "exporter.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "string.h"
#include "utility.h"
using namespace std;

const char template_text[] =
"image {\n"
"    resolution 900 600\n"
"    aa 0 2\n"
"    %  filter gaussian\n"
"}\n"

"trace-depths {\n"
"    diff 4\n"
"    refl 3\n"
"    refr 3\n"
"}\n"

"gi {\n"
"   type igi\n"
"   samples 64         % number of virtual photons per set\n"
"   sets 1             % number of sets (increase this to translate shadow boundaries into noise)\n"
"   b 0.00005          % bias - decrease this values until bright spots dissapear\n"
"   bias-samples 0     % set this >0 to make the algorithm unbiased\n"
"}\n"

"camera {\n"
"    type pinhole\n"
"    eye    50 -130 100\n"
"    target 50 55 20\n"
"    up     0 0 1\n"
"    fov    45\n"
"    aspect 1.5\n"
"}\n"

"shader {\n"
"   name Liquid\n"
"   type shiny\n"
"   diff { \"sRGB nonlinear\" 0.8 0.9 1.0 }\n"
"   refl 0.5\n"
"}\n"

"object {\n"
"    shader none\n"
"    type cornellbox\n"
"    corner0 0 0 0\n"
"    corner1  100  100 100\n"
"    left    0.80 0.25 0.25\n"
"    right   0.25 0.25 0.80\n"
"    top     0.70 0.70 0.70\n"
"    bottom  0.70 0.70 0.70\n"
"    back    0.70 0.70 0.70\n"
"    emit    15 15 15\n"
"    samples 32\n"
"}\n"

"light {\n"
"    type point\n"
"    color { \"sRGB nonlinear\" 1.000 1.000 1.000 }\n"
"    power 90000\n"
"    p 50 -200 100\n"
"}\n";



static void writeSunflow( int frame, vector<double> &vertices, vector<double> &normals, vector<int> &faces, FLOAT wall_thick ) {
	// Create Surface
	double s = 100.0*1.0/(1.0-2.0*wall_thick);
    bool fisrtTime = true;
    
    // If first time, write template
    if( fisrtTime ) {
        fisrtTime = false;
        system( "mkdir render/sunflow" );
        FILE *fp = fopen( "render/sunflow/template.sc", "w" );
        if( ! fp ) {
            dump( "Couldn't Write Sunflow Template File.\n");
            return;
        }
        fwrite( template_text, 1, strlen(template_text), fp );
        fclose(fp);
    }
    
	// Write into a new scene file
	char path[64];
	sprintf(path, "render/sunflow/%d_scene.sc", frame );
    
	FILE *fp = fopen( path, "w" );
    if( ! fp ) {
        dump( "Couldn't Write Sunflow File At %s\n", path );
        return;
    }

	// Write Object Header
    fprintf( fp, "include template.sc\n\n" );
	fprintf( fp, "object { \n shader Liquid \n type mesh \n name polySurfac \n %d %d \n", vertices.size()/3, faces.size()/3 );
	
	// Write Vertices and Normals
	for( int i=0; i<vertices.size(); i+=3 ) {
        FLOAT h = 0.001;
		fprintf( fp, "v %lf %lf %lf %lf %lf %lf 0 0\n", 
                s*fmin(1.0-h,fmax(h,(vertices[i+0]-wall_thick))), 
                s*fmin(1.0-h,fmax(h,(vertices[i+2]-wall_thick))), 
                s*fmin(1.0-h,fmax(h,(vertices[i+1]-wall_thick))), 
                normals[i], normals[i+2], normals[i+1] );
	}
	
	// Write Faces
	for( int i=0; i<faces.size(); i+=3 ) {
		fprintf( fp, "t %d %d %d\n", faces[i], faces[i+1], faces[i+2] );
	}
	
	fprintf( fp, "}\n" );
	fclose(fp);
}

void write_obj( int frame, vector<double> &vertices, vector<double> &normals, vector<int> &faces, FLOAT wall_thick) {
	
    double s = 1.0/(1.0-2.0*wall_thick);
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
		fprintf( obj_fp, "v %lf %lf %lf\n", s*(vertices[i]-wall_thick), s*(vertices[i+2]-wall_thick), s*(vertices[i+1]-wall_thick) );
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
					   vector<Object> &objects, std::vector<particle *> &particles, FLOAT wall_thick ) {
	write_obj( step, vertices, normals, faces, wall_thick );
    writeSunflow( step, vertices, normals, faces, wall_thick );
}
