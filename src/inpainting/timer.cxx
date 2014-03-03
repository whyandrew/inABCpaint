#include "../vxl_includes.h"

std::clock_t time_main = 0;
std::clock_t time_gradient = 0;
std::clock_t time_lookup = 0;
std::clock_t time_normal = 0;
std::clock_t time_conf = 0;
std::clock_t time_compute = 0;
unsigned int count_gradient = 0;
unsigned int count_normal = 0 ;
unsigned int count_conf = 0;
unsigned int count_lookup = 0;