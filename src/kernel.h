#include "common.h"

#ifndef _KERNEL_H_
#define _KERNEL_H_

namespace kernel {
    FLOAT smooth_kernel( FLOAT r2, FLOAT h );
    FLOAT sharp_kernel( FLOAT r2, FLOAT h );
}

#endif