#include <fstream>
#include<string>
#include<sstream>
#include "my_lib.h"
#include "lib.h"
#include <itpp/itbase.h>
using namespace itpp;
using namespace std;

string NumberToString(int pNumber) //not used anywhere. 
{
  ostringstream o;
  o<<pNumber;
  return o.str();
}

GF2mat append_vector(GF2mat G,bvec b){
  //append a row vector to the GF2mat G
  //used when saving error vectors for decoding
  G.set_size(G.rows()+1,G.cols(),true);
  G.set_row(G.rows()-1,b);
  return G;
}

GF2mat get_GF2mat(char * filename_prefix, char * filename_suffix){
  char filename[255];
  sprintf(filename,"%s%s",filename_prefix,filename_suffix);
  GF2mat E=MM_to_GF2mat(filename);
  return E;
}

GF2mat get_GF2mat(char * parent_folder, char * folder,  char *filename){
  char filename_full[255];
  sprintf(filename_full,"%s/%s/%s",parent_folder, folder, filename);
  GF2mat E=MM_to_GF2mat(filename_full);
  return E;
}

double get_error_density(GF2mat E){
  //fisrt row is zero
  //remove first row and return density of the remained submatrix
  double w;
  if(E.rows()==1){
    w=0;
  }else{
    w=E.get_submatrix(1,0,E.rows()-1,E.cols()-1).density();
  }
  return w;
}
