#ifndef MM_READ_H
#define MM_READ_H

#include <itpp/itbase.h>
//#include <string>
using namespace std;
using namespace itpp;

GF2mat MM_to_GF2mat(char *  file_name);
mat MM_to_mat(char *  file_name);

#endif
