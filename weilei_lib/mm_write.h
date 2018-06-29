#ifndef MM_WRITE_H
#define MM_WRITE_H

#include <itpp/itbase.h>
using namespace itpp;

//int GF2mat_to_MM(GF2mat G, char* file_name="mm_temp.dat");
int GF2mat_to_MM(GF2mat G, char* file_name);

//int mat_to_MM(mat G, char* file_name="mm_temp.dat");
int mat_to_MM(mat G, char* file_name);

#endif

