/*Weilei Zeng, 
March 18
from previous 20 ses of data for rand decoder with perm =5, accumulate them to get result with perm ={5..100..5}
 
check if the failure rate has exponential decay on # of permutation trials

March 2018,
make a z error-only decoder. G_x G_z seperate

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


double decode(char * filename_G, char * filename_E_input, char * filename_E_output1, char * filename_E_output2, char * filename_E_out, double p);

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

int main(int argc, char ** argv){
  char * filename_G=argv[1];//G_x
  char * filename_E_input=argv[2];//input z errors
  char * filename_E_output1=argv[3];//output z errors
  char * filename_E_output2=argv[4];//output z errors
  char * filename_E_output=argv[5];//output z errors
  double p=atof(argv[6]);//convert chars to integers
  p=p/100000.0;//1000.0;//because linux bash can only do intergers
  double failure_rate =
    decode(filename_G, filename_E_input,
	   filename_E_output1,filename_E_output2,
	   filename_E_output,p);
  return 0;
}


//double decode(char * filename_G, char * filename_S, char * filename_E, char * filename_E_out, double p){//p is the error probability for the depolarizing channel
double decode(char * filename_G, char * filename_E_input, char * filename_E_output1, char * filename_E_output2, char * filename_E_out, double p){//p is the error probability for the depolarizing channel
  //output the failure rate

  GF2mat G=MM_to_GF2mat(filename_G);//G_z
  int N=G.cols();///2;//number of qubits, size of the lattice
  G.set_size(G.rows()+1,G.cols(),true);//add one row, use G here because we will use diff=e+X later
 
  bvec e_t(N), e_d1(N),e_d2(N), e_d(N), diff_t(N);//e_t(2*N);//e_tilde=e_z//(e_z,e_x)
  //read error from file
  GF2mat E_input=MM_to_GF2mat(filename_E_input);
  
  int e_try=E_input.rows()-1;//number of total input errors
  //   e_try=(10<e_try) ? 10 :e_try; limit it to 10 to smaller
  GF2mat E_output(e_t,false),
    E_output1=MM_to_GF2mat(filename_E_output1),
    E_output2=MM_to_GF2mat(filename_E_output2),
    E_input_good(e_t,false),E_input_bad(e_t,false),
    E_output_good(e_t,false),E_output_bad(e_t,false);//false for row vectors.
  
  int e_bad=0;//count of bad errors
  for(int i1=0;i1<e_try;i1++){
      e_t=E_input.get_row(i1+1);//get input error
      e_d1=E_output1.get_row(i1+1);//get output error
      e_d2=E_output2.get_row(i1+1);//get output error
      
      if (weight(e_d1)>weight(e_d2) ){
	  e_d=e_d2;
	}else{
	  e_d=e_d1;
	}
	//check the error and syndrome    
      bvec diff_t=e_d+e_t;//same format, (e_z,e_x)

      E_output = append_vector(E_output, e_d);

      G.set_row(G.rows()-1,diff_t);//check if the diff belong to gauge group
      //cout<<G.rows()<<", rank of new G = "<<G.row_rank()<<endl;
      if(G.rows()==G.row_rank()){//not belond to gauge -> belongs to logical group -> bad error
	e_bad++;
	//	cout<<"BAD* ";
	//cout<<"weight(e_t) = "<<weight(e_t)<<", wmin = "<<wmin<<", weight(diff_t) = "<<weight(diff_t)<<endl;
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

    //save error output

    //    char filename_o[255];
    //    sprintf(filename_o, "%s_output",filename_E_out);
    GF2mat_to_MM(E_output,filename_E_out);
   
    
    char filename_ig[255];
    sprintf(filename_ig, "%s_input_good",filename_E_out);
    GF2mat_to_MM(E_input_good,filename_ig);
    char filename_ib[255];
    sprintf(filename_ib, "%s_input_bad",filename_E_out);
    GF2mat_to_MM(E_input_bad,filename_ib);
    char filename_og[255];
    sprintf(filename_og, "%s_output_good",filename_E_out);
    GF2mat_to_MM(E_output_good,filename_og);
    char filename_ob[255];
    sprintf(filename_ob, "%s_output_bad",filename_E_out);
    GF2mat_to_MM(E_output_bad,filename_ob);
    
    //print result
    cout<<"# of bonds/qubits N = "<<N
	<<", # of total input error e_try= "<<e_try
	<<", p = "<<p
	<<", failure_rate = "<<failure_rate
	<<endl;
    
    return failure_rate;
}
