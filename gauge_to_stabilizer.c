/*Weilei Zeng, Nov,2017
For a subsystem code, from gauge generating matrix G, calculate S and save it in .mm file.
 */

#include <stdio.h>
//#include <stdlib.h>
#include<itpp/itbase.h>
#include <itpp/itcomm.h>
#include "my_lib.h"
using namespace itpp;
using namespace std;

//int decoding();
int G_to_S(char * filename_G, char * filename_S);

GF2mat get_tilde(GF2mat G){//do I need to used &  ?
  GF2mat G_t=G.get_submatrix(0,G.cols()/2,G.rows()-1,G.cols()-1).concatenate_horizontal(  G.get_submatrix(0,0,G.rows()-1,G.cols()/2-1) );
  //get G_tilde, G=(G_x,G_z), G_tilde=G_t=(G_z,G_x);
  return G_t;
}
int weight(bvec &b){//we can use BPSK::count(zero_b(N),b) to count the weight
  int length=b.length();
  int weight=0;
  for (int i=0;i<length;i++){
    if(b.get(i)){
      weight++;
    }
  }
  return weight;
}
//int decode_case1();
//int decode_case2();
int main(int argc, char** argv){
  char * filename_G=argv[1];
  char * filename_S=argv[2];
  G_to_S(filename_G,filename_S);
  return 0;
}

/*
int decode_case1(){//case
  //decode the corresponding file in data/
  //first run n=(5,23,2) finished
  //second run n=(23,33,2) running
  //thrid run n=(33,43,2)
  for (int n=33;n<43;n+=2){
    for(double d=0.5;d<0.91;d+=0.1){
      for(int i=0;i<10;i++){
	char filename_A[66],filename_G[66];
	sprintf(filename_G,"data/bravyi_%s_size_%d_density_%.2f_v_%d.dat","G",n,d,i);
	char filename_H[66],filename_Q[66];
	sprintf(filename_H,"data/bravyi_%s_size_%d_density_%.2f_v_%d.dat","H",n,d,i);
	sprintf(filename_Q,"data/bravyi_%s_size_%d_density_%.2f_v_%d.dat","Q",n,d,i);
	//	cout<<filename_G<<endl<<filename_H<<endl<<filename_Q<<endl;
	//	decode(filename_G,filename_H,filename_Q);
      }
    }
  }
  return 0;
}
*/

int G_to_S(char * filename_G, char * filename_S){

  //set up parameters
  //  double p_initial=0.010, p_final=0.050, p_increment=0.002;
  // int intervals=1+floor(0.1+(p_final-p_initial)/p_increment );
  // int e_try=100, perm_try=5;
  //mat mat_data(2,intervals);//save p and failure rate;
  
  //calculate stabilzier S from gauge G. run over all errors and determine its good or bad.

  //cout<<filename_G<<endl;
  //  GF2mat G=mm_read(filename_G);
  GF2mat G=MM_to_GF2mat(filename_G);
  int N=G.cols()/2;//number of qubits, size of the lattice
  //cout<<"finish reading G"<<endl;
  //  cout<<"G="<<endl<<G<<endl;

  //use Gaussian elimination to find S, because S is an abelian subgroup of G
  GF2mat G_t=get_tilde(G);
  //get G_tilde, G=(G_x,G_z), G_tilde=G_t=(G_z,G_x);
  //cout<<"G_t="<<endl<<G_t<<endl;
  GF2mat ss=G_t*(G.transpose());
  GF2mat T,U;
  ivec P;//used for all the gaussian elimination in this file, T_fact(T,U,P)
  ss.T_fact(T,U,P);
  //cout<<"T="<<T<<endl<<"U="<<U<<endl<<"P="<<P<<endl;
  //sometimes I get a full rank of U here, means there is no stabilizer for this G
  int rank_U=U.row_rank();
  if(rank_U==U.rows()){
    cout<<"No stabilzier for this G. Abort the calculation. Filename: "<<filename_G<<endl;
    return 0;
  }
  GF2mat S=(T*G).get_submatrix(rank_U,0,G.rows()-1,G.cols()-1);
  //cout<<"get S and H"<<endl;
  cout<<"S="<<endl<<S<<endl;
  GF2mat_to_MM(S,filename_S);

  //cout<<"test commutation: G_t*(S.transpoze()) is "<<G_t*(S.transpose())<<endl;

  return 0;
}
