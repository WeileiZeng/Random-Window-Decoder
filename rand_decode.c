/*Weilei Zeng, test file on 3/13/2018
March 2018,
make a z error-only decoder. G has only X operators.

Feb,2018
Compare to rand_decode.c, this file is design for bp_decoding3.c. The output errors from bp decoding is saved and used as the input error for the rand decoder here.

old notes: Nov,2017
For a subsystem code, given G and S, in a depolorizing channel with error probability p, generating errors can decode it. The output is the ratio of decoding failure.
Decoding method: Random Window Decoder. The detail was explained in the PDF file "weilei reasearch note".
This small program is aim to seperate the derivation of S and the decoding, which works in parallel in linux, so that we can take advantage of the mulitple CPUs.
speed: the time to do one permutation on one error for 25 qubit bacon shor code is 22/7=3 seconds.
10,000 permutation means 30,000 sec=10 CPU hours

 */

#include <stdio.h>
#include<itpp/itbase.h>
#include <itpp/itcomm.h>
#include "my_lib.h"
using namespace itpp;
using namespace std;

//int decoding();
double decode(char * filename_G, char * filename_S, char * filename_E, double p);
//double decode(char * filename_G, char * filename_S, double p);
/*donot need in this z error only decoder. for me it looks like a classical decoder
GF2mat get_tilde(GF2mat G){//do I need to used &  ?
  //this function is not used in this program
  GF2mat G_t=G.get_submatrix(0,G.cols()/2,G.rows()-1,G.cols()-1).concatenate_horizontal(  G.get_submatrix(0,0,G.rows()-1,G.cols()/2-1) );
  //get G_tilde, G=(G_x,G_z), G_tilde=G_t=(G_z,G_x);
  return G_t;
}
bvec get_tilde(bvec e){//do I need to used &  ?
 //get e_tilde, e=(e_x,e_z), e_tilde=e_t=(e_z,e_x);or the other way around
  int len=e.length();
  bvec e_t(len);
  e_t.set_subvector(0,e.get(len/2,-1) );
  e_t.set_subvector(len/2,e.get(0,len/2-1) );
  return e_t;
}
*/

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
/*
int save_result(double p, double failure_rate,char * filename){
  mat mat_data(2,1);
  mat_data.set(0,0,p);
  mat_data.set(1,0,failure_rate);
  mat_to_MM(mat_data,filename);
  return 0;
}
*/
int main(int argc, char ** argv){
  char * filename_G=argv[1];
  char * filename_S=argv[1];//S and G are the same for stabilizer code, like toric code //leave it here for future compatibility
  char * filename_E=argv[2];//input errors
  //char * filename_result=argv[3];
  double p=atof(argv[3]);//convert chars to integers
  p=p/100000.0;//1000.0;//because linux bash can only do intergers
    //char filename_result_p[255];
  //sprintf( filename_result_p,"%s%.5f", filename_result,p);//division 0.00001
  //   cout<<"f r ="<<filename_result_p<<endl;
  double failure_rate = decode(filename_G,filename_S,filename_E,p);
  //  save_result(p,failure_rate,filename_result_p);
  return 0;
}


double decode(char * filename_G, char * filename_S, char * filename_E, double p){//p is the error probability for the depolarizing channel
  //output the failure rate
  RNG_randomize();//get randome seed 
  //#set up parameters#
  // int e_try=100;//number of random errors generated//get from E
  int perm_try=5;//number of trails of random window / permutation; It seems that the minmum weight can be reached mostly within 5 times for 25 qubit bacon shor code.//5 or 7 for toric 5x5 has no difference. 15 make it slightly better
  GF2mat G=MM_to_GF2mat(filename_G);
  int N=G.cols();///2;//number of qubits, size of the lattice
  //  cout<<"G="<<endl<<G<<endl;
  GF2mat T,U;
  ivec P;//used for all the gaussian elimination in this file, T_fact(T,U,P)
  GF2mat S=MM_to_GF2mat(filename_S);
  //cout<<"S="<<endl<<S<<endl;

  //cout<<"test commutation: G_t*(S.transpoze()) is "<<G_t*(S.transpose())<<endl;
  //define H=(S_x|S_z|s)
  //find Q, the dual of H, which is e+G+L, HQ^T=0
  
  GF2mat H=S;    
  H.set_size(S.rows(),S.cols()+1,true);//H=(H_x,s)//H=(H_x,H_z,s)
  G.set_size(G.rows()+1,G.cols(),true);//add one row, use G here because we will use diff=e+X later
 
  bvec e_t=zeros_b(N);//e_t(2*N);//e_tilde=e_z//(e_z,e_x)
  //read error from file
  GF2mat E_input=MM_to_GF2mat(filename_E);//here the input is the output of the nonconvergent cases after BP decoding. The first row of this vector is an extra zero vector.
  int e_try=E_input.rows()-1;//number of total input errors
  e_try=(10<e_try) ? 100 :e_try;// limit it to 10 to smaller
  GF2mat E_input_good(e_t,false),E_input_bad(e_t,false),E_output_good(e_t,false),E_output_bad(e_t,false);//false for row vectors.
  
  //  ivec perm_half(S.cols()/2), perm(S.cols()+1);//permutation vectors
  ivec perm(S.cols()+1);
  perm.set(S.cols(),S.cols());//the last col for syndrome, is fixed

  int e_bad=0;//count of bad errors
  for(int i1=0;i1<e_try;i1++){
    //      for (int i2=0;i2<2*N;i2++){//setup random error with error rate p
    //	e_t.set(i2,(randu()-p<0)? 1:0); 
    // }

      e_t=E_input.get_row(i1+1);//get input error
      bvec s=S*e_t;//syndrome s=s_x;//(s_x,s_z);
      H.set_col(H.cols()-1,s);//add syndrome
      //cout<<"parity check matrix H=(S_x,S_z,s). The last column is the syndrome. "<<H<<endl;
    
      int wmin=e_t.length();//minimum weight
      bvec e_d=zeros_b(e_t.length() );//error detected with minimum weight wmin //e_d =(e_z,e_x)
      // int perm_try=20;//number of trials to decode the error
      for (int i3=0;i3<perm_try;i3++){//find the error with minimum weight
	//set up random permutation vector for H
	perm.set_subvector(0,sort_index(randu( S.cols() )) );
	//cout<<perm<<endl;
	/*
	perm_half=sort_index( randu( S.cols()/2 ) );
	perm.set_subvector(0,perm_half);
	perm.set_subvector(S.cols()/2,perm_half+S.cols()/2);
	//use the same permutation for H_x and H_z
	*/

	H.permute_cols( perm, false);
	//better to check the commutation again for H
	H.transpose().T_fact(T,U,P);//can I use same T,U,P matrix/vec here despite different size of them: Yes, any size and any value won't affect the result
	//The rows of this matrix are respectively: all errors with non-zero syndrome; all errors with zero syndrome (gauge errors and logical errors); error with the given syndrome. This maybe useful.

	//T should include all the orther inequivelant errors. The resize make it cleaner but may lose some meaningful info
	GF2mat Q=T.get_submatrix(H.rows(),0,H.cols()-1,H.cols()-1);//does permutation matter?  It doesn't matter cause we only look at the zero rows in the H_z.transpose() after gauss. I don't need the nonzero rows to be triangle or not
	//In fact, this Q=(Q_z,Q_x,1) should be Q_tilde. But it is up to the definition, and doesn't affect the permutation
	//cout<<(H*Q.transpose()).density()<<endl;
	//cout<<"rank, row = "<<Q.row_rank()<<", "<<Q.rows()<<endl;
	//cout<<(H*Q.transpose())<<endl;
	Q.permute_cols(perm,true);//permute back
	/*check permutation relation of Q,H,G
	cout<<"G*Q^T, "<<
	  (G.get_submatrix(0,0,G.rows()-2,G.cols()-1)*
	  Q.get_submatrix(0,0,Q.rows()-2,Q.cols()-2).transpose()
	   ).density()
	   <<endl;
	cout<<"rank of G: "<<G.get_submatrix(0,0,G.rows()-2,G.cols()-1).T_fact(T,U,P)<<endl;
	cout<<"rank of Q: "<<Q.get_submatrix(0,0,Q.rows()-2,Q.cols()-2).transpose().T_fact(T,U,P)<<endl;
	*/
	
	//	cout<<"Q, "<<Q<<endl;
	//cout<<"Q: The last column means if the error match zero syndrome or match the given syndrome in H. Hence only the last row is the error we want."<<endl;

	//get error and check
	//Q.T_fact(T,U,P);//for test, not needed//seems no improvement
	bvec X_t=Q.get_row(Q.rows()-1);//the error detected, X_t=(X_z,X_x,1)

	//get the row with minimum weight. add X_t to make sure the last element is one
	int ww=weight(X_t);
	for (int a = 0;a<Q.rows()-1;a++){
	  if (ww> weight(X_t+Q.get_row(a)) ){
	    ww= weight(X_t+Q.get_row(a));
	    X_t=X_t+Q.get_row(a);
	    if (Q.get_row(a)(X_t.length()-1) ==1) cout<<"Ding"<<endl;
	    if (X_t(X_t.length()-1)==0)cout<<"damn it!"<<endl;
	  }
	}
	X_t.set_size(X_t.size()-1,true);//remove the last 1.  X_t=(X_z,X_x)
	int w=weight(X_t);
	//cout<<"weight = "<<w<<endl;
	if(wmin>w){//find another error with smaller weight, update it
	  wmin=w;
	  e_d=X_t;
	}
	H.permute_cols(perm, true);//permute back for another run
      }
  
	//check the error and syndrome    
      bvec diff_t=e_d+e_t;//same format, (e_z,e_x)
      cout<<"weight(e_t) = "<<weight(e_t)<<", wmin = "<<wmin<<", weight(diff_t) = "<<weight(diff_t)<<endl;
      // bvec diff=get_tilde(diff_t);
      //      G.set_row(G.rows()-1,diff);//check if the diff belong to gauge group
      G.set_row(G.rows()-1,diff_t);//check if the diff belong to gauge group
      //cout<<G.rows()<<", rank of new G = "<<G.row_rank()<<endl;
      if(G.rows()==G.row_rank()){//not belond to gauge -> belongs to logical group -> bad error
	e_bad++;
	//cout<<"BAD*";
	E_input_bad = append_vector(E_input_bad,e_t);
	E_output_bad = append_vector(E_output_bad,e_d);
	
      }else{
	//cout<<"__good__";
	E_input_good = append_vector(E_input_good,e_t);
	E_output_good = append_vector(E_output_good,e_d);
	
      }
      //cout<<S*diff_t<<endl;//check if get the same syndrome
  }
    double failure_rate=1.0*e_bad/e_try;

    //cout<<"E_output_bad "<<E_output_bad<<endl;
    //cout<<"p="<<p<<", failure_rate="<<failure_rate<<endl;
    //cout<<"number of input errors:"<<e_try<<endl;
    //cout<<"number of bad errors:"<<e_bad<<endl;

    /*
    //save error output
    char filename_ig[255];
    sprintf(filename_ig, "%s_input_good",filename_E);
    GF2mat_to_MM(E_input_good,filename_ig);
    char filename_ib[255];
    sprintf(filename_ib, "%s_input_bad",filename_E);
    GF2mat_to_MM(E_input_bad,filename_ib);
    char filename_og[255];
    sprintf(filename_og, "%s_output_good",filename_E);
    GF2mat_to_MM(E_output_good,filename_og);
    char filename_ob[255];
    sprintf(filename_ob, "%s_output_bad",filename_E);
    GF2mat_to_MM(E_output_bad,filename_ob);
    */
    
    //print result
    cout<<"# of bonds/qubits N = "<<N
	<<", # of total input error e_try= "<<e_try
	<<", p = "<<p
	<<", failure_rate = "<<failure_rate
	<<endl;
    
    return failure_rate;
}
