#include "flip3D.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char * argv[]) {
	int frame = 0;
	if( argc == 2  ) {
		sscanf( argv[1], "%d", &frame );
		printf( "frame = %d\n", frame );
	}
	flip3D::init(frame);
	while(1) flip3D::simulateStep();
	return 0;
}
