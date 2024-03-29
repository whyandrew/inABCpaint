//  DO NOT MODIFY THIS FILE!!
//
//
// Include files from the VXL library

// Files that replace standard C++ libraries
#include <vcl_cstdlib.h>
#include <vcl_cstdio.h>
#include <vcl_cassert.h>
#include <vcl_string.h>
#include <vcl_iostream.h>
#include <vcl_queue.h>

// The libraries below are required to manipulate VXL images
#include<core/vil/vil_load.h>
#include<core/vil/vil_save.h>
#include<core/vil/vil_image_view.h>
#include<core/vil/vil_rgb.h> 
#include<core/vil/vil_convert.h>
#include<core/vil/vil_math.h>
#include<core/vil/vil_crop.h>
#include<core/vil/vil_resample_bilin.h>
#include<core/vil/vil_clamp.h>
#include"vil_trace_8con_boundary.h"

#include<core/vul/vul_arg.h>

// VXL Libraries used for numerical computations
//#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_double_2x2.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_2.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_rational.h>

#include <ctime>

#ifndef TIMER_H
#define TIMER_H

#define TIMER_START     std::clock_t start_time;	start_time = std::clock()
#define TIMER_ELLAPSED  std::clock() - start_time

    extern std::clock_t time_main;
    extern std::clock_t time_gradient;
    extern std::clock_t time_lookup;
    extern std::clock_t time_normal;
    extern std::clock_t time_conf;
    extern unsigned int count_gradient;
    extern unsigned int count_normal;
    extern unsigned int count_conf;
    extern unsigned int count_lookup;

#endif