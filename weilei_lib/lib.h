#ifndef LIB_H
#define LIB_H
#include <string>
#include <itpp/itbase.h>
using namespace itpp;
using namespace std;


string NumberToString(int pNumber);
GF2mat append_vector(GF2mat G,bvec b);
GF2mat get_GF2mat(char * filename_prefix, char * filename_suffix);
GF2mat get_GF2mat(char * parent_folder, char * folder, char * filename);
double get_error_density(GF2mat E);
#endif
