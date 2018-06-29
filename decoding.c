/*Weilei Zeng, Nov,2017
For a subsystem code, from gauge generating matrix G, calculate S. Then, in a depolorizing channel with error probability p, generating errors and then decode it. The output is the ratio of decoding failure.
Decoding method: Random Window Decoder. The detail was explained in the PDF file "weilei reasearch note".
This program do all the calculation together. But it is easier to seperate the derivation of S and the decoding, which then can work in parallel in linux. Hence, we can take advantage of the mulitple CPUs.
 */

#include <stdio.h>
//#include <stdlib.h>
#include<itpp/itbase.h>
#include <itpp/itcomm.h>
#include "my_lib.h"
using namespace itpp;
using namespace std;

int decoding();
int decode(char * filename_G, char * filename_H, char* filename_Q);

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
int decode_case1();
int decode_case2();
int main(){
  decode_case2();
  //  decoding();
  //  decode("G_BS_size5.dat","tempH.dat","tempQ.dat");
  return 0;
}
int decode_case2(){
  decode("data/bacon/BaconShor_size_9_G.mm","mat_data_9.mm","tempQ.dat");
  //  decode("data/bacon/BaconShor_size_25_G.dat","mat_data_25.dat","tempQ.dat");
  //decode("data/bacon/BaconShor_size_51_G.dat","mat_data_51.dat","tempQ.dat");

  return 0;
}


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
	decode(filename_G,filename_H,filename_Q);
      }
    }
  }
  return 0;
}


int decode(char * filename_G, char * filename_H, char* filename_Q){

  //set up parameters
  double p_initial=0.010, p_final=0.050, p_increment=0.002;
  int intervals=1+floor(0.1+(p_final-p_initial)/p_increment );
  int e_try=100, perm_try=5;
  mat mat_data(2,intervals);//save p and failure rate;
  
  //calculate stabilzier S from gauge G. run over all errors and determine its good or bad.

  //cout<<filename_G<<endl;
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
  //  cout<<"S="<<endl<<S<<endl;
  //GF2mat_to_MM(S,filename_S);

  //cout<<"test commutation: G_t*(S.transpoze()) is "<<G_t*(S.transpose())<<endl;
  //define H=(Q|s)
  //find Q, the dual of H, which is e+G+L, HQ^T=0
   
  GF2mat H=S;
  RNG_randomize();//get randome seed 
  //  cout<<"permutation -->perm_half:"<<perm_half<<endl;
  //cout<<"permutation -->perm:"<<perm<<endl;
  //use the same permutation for H_x and H_z
    
  H.set_size(S.rows(),S.cols()+1,true);//H=(H_x,H_z,s)

  //cout<<"rank of G"<<G.row_rank()<<"  rows of G"<<G.rows()<<endl;
  G_t.set_size(G.rows()+1,G.cols(),true);//add one row, use G_t here because we will use e_t and X_t later
  //  G.set_row(G.rows()-1,e);//initialize it for nonsense, guess it should have been zero there
 
  bvec e_t(2*N);//e_tilde=(e_z,e_x)//=randb(2*N);//random error p=0.5
  

  ivec perm_half(S.cols()/2);//=sort_index( randu( S.cols()/2 ) );
  ivec perm(S.cols()+1);
  perm.set(S.cols(),S.cols());//the last col for syndrome, is fixed

  double p=p_initial-p_increment;//error probability
  //  int intervals=25;

  for(int ip=0;ip<intervals;ip++){
    p+=p_increment;
    int e_bad=0;//number of bad errors
    //  int e_try=30;//number of trials for different errors  
    for(int i1=0;i1<e_try;i1++){
      for (int i2=0;i2<2*N;i2++){//setup random error
	e_t.set(i2,(randu()-p<0)? 1:0);
      }

      //      cout<<endl<<"e_t  ="<<e_t<<endl;
      //  bvec s=randb(S.rows());//="0 1 0 1";//sample syndrome, the length should match the number of stabilziers
      bvec s=S*e_t;//syndrome s=(s_x,s_z);
    
      //cout<<"The syndrome is "<<s<<endl;
      H.set_col(H.cols()-1,s);//add syndrome
      //cout<<"parity check matrix H=(S_x,S_z,s). The last column is the syndrome. "<<H<<endl;
    
      int wmin=e_t.length();//minimum weight
      bvec e_d=e_t;//error detected with minimum weight wmin
      // int perm_try=20;//number of trials to decode the error
      for (int i3=0;i3<perm_try;i3++){//find the error with minimum weight
	//set up random permutation vector for H
	perm_half=sort_index( randu( S.cols()/2 ) );
	perm.set_subvector(0,perm_half);
	perm.set_subvector(S.cols()/2,perm_half+S.cols()/2);
  
	H.permute_cols( perm, false);
	//  cout<<"H"<<H<<endl;
	//better to check the commutation again for H
  
	H.transpose().T_fact(T,U,P);//can I use same matrix here despite different size of them: Yes, any size and any value won't affect the result
	//The rows of this matrix are respectively: all errors with non-zero syndrome; all errors with zero syndrome (gauge errors and logical errors); error with the given syndrome. This maybe useful.
  
	//cout<<"H after gauss"<<endl<<T*(H_z.transpose())<<endl;
	// cout<<"T"<<endl<<T<<endl;

	//T should include all the orther inequivelant errors. The resize make it cleaner but may lose some meaningful info
	GF2mat Q=T.get_submatrix(H.rows(),0,H.cols()-1,H.cols()-1);//does permutation matter?  It doesn't matter cause we only look at the zero rows in the H_z.transpose() after gauss. I don't need the nonzero rows to be triangle or not
	//In fact, this Q=(Q_z,Q_x,1) should be Q_tilde. But it is up to the definition, and doesn't affect the permutation
	//      cout<<"Q"<<Q<<endl;
	Q.permute_cols(perm,true);//permute back

	//  cout<<"permute back Q=(Q_z,Q_x,1)"<<endl<<Q<<endl;
	//  cout<<"get Q"<<endl;
	//cout<<"Q: The last column means if the error match zero syndrome or match the given syndrome in H. Hence only the last row is the error we want."<<endl;

	//      GF2mat_to_MM(Q,filename_Q);

	//get error and check
	bvec X_t=Q.get_row(Q.rows()-1);//the error detected
	X_t.set_size(X_t.size()-1,true);//remove the last 1
	int w=weight(X_t);
	if(wmin>w){//find another error with smaller weight, update it
	  wmin=w;
	  e_d=X_t;
	}
	//check the error and syndrome
  
	H.permute_cols(perm, true);//permute back for another run
	//      GF2mat_to_MM(H,filename_H);
      }
    
      bvec diff_t=e_d+e_t;
      //    cout<<"e_d  ="<<e_d<<endl<<"diff ="<<diff_t<<endl;
      G_t.set_row(G_t.rows()-1,diff_t);//check if the diff belong to gauge group
      //      cout<<"rows "<<G.rows()<<"  rank"<<G.row_rank()<<endl;
      if(G_t.rows()==G_t.row_rank()){
	//cout<<"*BAD*";
	e_bad++;
      }else{
	//cout<<"good_";
      }
    
    
    }
    double failure_rate=1.0*e_bad/e_try;
    mat_data.set(0,ip,p);
    mat_data.set(1,ip,failure_rate);
    cout<<"p="<<p<<", failure_rate="<<failure_rate<<endl;
  }
  mat_to_MM(mat_data,filename_H);//"mat.dat");
    /*
  cout<<"test commutation: H*Q   "<<
    H*Q.transpose()
    //     H.get_submatrix(0,0,H.rows()-1,H.cols()-2)     * ( ( Q.get_submatrix( 0,0, Q.rows()-1, Q.cols()-2  )  ).transpose()  )
    //        get_tilde( H.get_submatrix(0,0,H.rows()-1,H.cols()-2)   )  * ( ( Q.get_submatrix( 0,0, Q.rows()-1, Q.cols()-2  )  ).transpose()  )
    <<endl;

  cout<<"Number of qubits             ----->"<<G.cols()/2<<endl
      <<"Number of gauge generator    ----->"<<G.rows()<<endl
      <<"Number of Stabilzier generators--->"<<S.rows()<<endl
      <<"Number of logical operators  ----->"<<Q.rows()-1-G.rows()<<endl;
  */  
  return 0;
}
