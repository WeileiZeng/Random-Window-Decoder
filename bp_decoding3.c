//In order to change to z error only, there is no change needed for this code. Just use the new generating matrix. 2/10/2018
//Compared to bp_decoding2.c, in this file,I save all data for errors that converge or not converge in bp decoding.
//In this file, for toric code, we use stabilizer S=(S_x,S_z).
//errors are in the format e=(e_z,e_x). (bistin, recbits, bitsout) //not checked excatly, based on the assumption that the BP decoding in itpp doesnot know the difference between X and Z errors and operators.
//hello world

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <sstream>
#include <stdio.h>
#include "my_lib.h"
using namespace std;
using namespace itpp;

/*
//not used in this file
int print_toric_error(bvec error_bits){//print error for the toric code
  //in stabilizer generating matrix, we have vertex for X and plaquette for Z
  //in this print, we print X error first and then Z errors to the right
  //use plaquette to check X error and vertex to check Z errors
  int N=error_bits.length();
  int n=(int) sqrt(N/4);
  cout<<"n="<<n<<endl;
    //      <<"N="<<N<<endl;
  for (int i=0;i<n;i++){
    
    //cout<<error_bits.get(0,N-1)<<endl;//start from 0 to N-1
    cout<<"h "<<error_bits.get(n*i,n*i+n-1)<<"   "
      	<<""<<error_bits.get(n*i+n*n*2,n*n*2+n*i+n-1)<<"   "<<endl
	<<"v"<<error_bits.get(n*i+n*n,  n*n+n*i+n-1  )<<"   "
	<<""<<error_bits.get(n*i+n*n*3,n*n*3+n*i+n-1)<<"   "
	<<endl;
  }
  
  return 0;

}
*/

LDPC_Code MM_to_LDPC_Code(char * filename){
  //convert GF2mat saved in .mm file to LDPC_Code
  GF2mat G=MM_to_GF2mat(filename);
  GF2mat_sparse Gs=G.sparsify();
  GF2mat_sparse_alist Gsa;
  Gsa.from_sparse(Gs);
  LDPC_Parity H(Gsa);
  LDPC_Code C(&H);
  return C;
}

int save_result(double p, double converge_rate, char * filename){
  mat mat_data(2,1);
  mat_data.set(0,0,p);
  mat_data.set(1,0,converge_rate);
  mat_to_MM(mat_data,filename);//dosen't work here and not sure why. error recorded in mm_write.c
  return 0;
}


// Read the code from files and do BP decoding
//input:source file for stabilzier matrix; error propability p ;

int main(int argc, char **argv){
  //parameter setup
  int cycles=10000;//number of cycles: fro toric code, 10000 give reletively clear result
  char * filename_G=argv[1];
  char * filename_result=argv[2];//prefix for the file
  double p=atof(argv[3]);
  p=p/100000.0;//previous use 1000 division. Now use 100,000 division cause the thershold for toric codes seems to be around 0.1%.
  char filename_result_p[255];
  sprintf( filename_result_p,"%s%.5f",filename_result,p);//append p in the file name
  
  Real_Timer timer;

  LDPC_Code C=MM_to_LDPC_Code(filename_G);  //convert GF2mat saved in .mm file to LDPC_Code

  //LDPC_Code C(filename_G);//load code if saved in .it file with LDPC_Code format
  //vec pvals = "0.01:0.01:0.1"; p.get(pvals,"pvals");// start:increment:end=inclusive

  //bp decoding set up
  // C.set_exit_conditions(2500);//high perperformance
  C.set_exit_conditions(50,true,true);  // 50 iterations,check syndrome always
  C.set_llrcalc(LLR_calc_unit(12,0,7)); // fast mode; see "llr.h" 
  //  cout << C << endl;
  int N = C.get_nvar(); // number of bits per codeword; N=2 x n x n is the size of the toric code.
  
  bvec bitsin = zeros_b(N);//original zero vector
  GF2mat E_input_converge(bitsin, false),E_input_nonconverge(bitsin, false),E_output_converge(bitsin,false),E_output_nonconverge(bitsin, false);//false for row vector;//These two GF2mat saves the input errors and their output after bp_decoding. They are clasified into convergeds cases and non converged cases.
  RNG_randomize();
  timer.tic(); 

  double pp=p;//=0.05;//test error probability //=pvals(j);   //  double pp=pow(10.0,-pvals(j));
  //pp=0.01  ans~1;pp=0.03  ans=-2500
  BSC bsc(pp); // initialize BSC, error channel with probability pp
  int ans=0;
      
  // bvec rec_bits = bitsin;   // input vector with manually-input errors

  int counts_converge=0;
  for (int i=0;i<cycles;i++){
    bvec rec_bits = bsc(bitsin);    // input vector with errors from bsc chanel
    //cout<<"rec_bits="<<rec_bits<<endl;
    vec s(N), s0(N); // input LLR version 
    QLLRvec llr;
    for(int i1=0;i1<N;i1++)
      s(i1)=s0(i1)=(rec_bits(i1)==1)?-log(1.0/pp-1.0):log(1.0/pp-1.0);//convert to input format for bp_decode  
    //cout<<"pp="<<pp<<"\t log(1.0/pp-1.0)="<<log(1.0/pp-1.0)<<endl;
    //cout<<"s="<<s<<endl;
  
    QLLRvec llr_input=C.get_llrcalc().to_qllr(s);
    //cout<<"llr_input="<<llr_input<<endl;
    ans=C.bp_decode(llr_input, llr);
    //cout<<"ans="<<ans<<endl;
    bvec bitsout = llr < 0;//output
    if(ans>=0){
      //  cout<<"llr="<<llr<<endl;
      //cout<<"bitsout="<<bitsout<<endl;
      //cout<<"diff = "<<bitsin+bitsout<<endl;
      //    cout<<"diff = "<<rec_bits+bitsout<<endl;
      //check min weight of bitsout?
      E_input_converge = append_vector(E_input_converge, rec_bits);
      E_output_converge = append_vector(E_output_converge, bitsout);//maybe all zero, but no harm to check            
      counts_converge++;
    }else{
      //            cout<<"not converge when i="<<i<<endl
      //  <<"rec_bits= "<<rec_bits<<endl
      //  <<"bitsout = "<<bitsout<<endl;
      //  cout<<"NC, \# of errors in rec_bits = "<<BERC::count_errors(bitsin,rec_bits);
      // cout<<", \# of errors in bitsout  = "<<BERC::count_errors(bitsin,bitsout)<<endl;
      //      print_toric_error(rec_bits);
      //print_toric_error(bitsout);

      E_input_nonconverge = append_vector(E_input_nonconverge,rec_bits);
      E_output_nonconverge = append_vector(E_output_nonconverge,bitsout);//this would be the input error for the random window decoder.
      // break;//break here to just print the first non converged error
    }
  }
  //cout<<"counts_converge="<<counts_converge<<endl;
  double rate_converge=counts_converge*1.0/cycles;
  cout<<"N = "<<N<<" =  2 x "<< N/2;
  cout<<", p = "<<p;
  cout<<", Converge rate ="<<rate_converge<<endl;
  //  save_result(p,rate_converge,filename_result_p);//no need to save this. can get it by counting the size of the matrix when doing gnuplot
  
  //save errors to files
  char filename_result_ic[255];
  sprintf( filename_result_ic,"%s%.5f_input_converge",filename_result,p);//append p to the file name
  //cout<<E_input_converge.rows()<<endl;
  GF2mat_to_MM(E_input_converge,filename_result_ic);  
    
  char filename_result_in[255];
  sprintf( filename_result_in,"%s%.5f_input_nonconverge",filename_result,p);//append p to the file name
  GF2mat_to_MM(E_input_nonconverge,filename_result_in);  

  char filename_result_oc[255];
  sprintf( filename_result_oc,"%s%.5f_output_converge",filename_result,p);//append p to the file name
  GF2mat_to_MM(E_output_converge,filename_result_oc);  

  char filename_result_on[255];
  sprintf( filename_result_on,"%s%.5f_output_nonconverge",filename_result,p);//append p to the file name
  GF2mat_to_MM(E_output_nonconverge,filename_result_on);  
  
  //  cout<<E_input_nonconverge<<endl;
  // cout<<E_input_converge<<endl;
  timer.toc_print();
  return 0;
}
   
