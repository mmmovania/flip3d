#include "kernel.h"
#include <math.h>

FLOAT kernel::smooth_kernel( FLOAT r2, FLOAT h ) {
    return fmax( 1.0-r2/(h*h), 0.0 );
}

FLOAT kernel::sharp_kernel( FLOAT r2, FLOAT h ) {
    return fmax( h*h/fmax(r2,1.0e-5) - 1.0, 0.0 );
}