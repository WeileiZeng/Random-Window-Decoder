#include "my_lib.h"
#include <itpp/itcomm.h>
#include <stdio.h>
using namespace std;
using namespace itpp;

int bp_converge2(){//result collection for bp decoding //100,000 division for p
  //  cout<<"another thing"<<endl;
  //size 5,9,25,35
  char filename_prefix[]="data/toric/bp_converge2/toric_S_size_9.mm";
  mat mat_data;
  for (double p=0.00010;p<0.00200;p+=0.00010){
    char filename_p[255];
    sprintf(filename_p,"%s%.5f",filename_prefix,p);
    cout<<filename_p<<endl;
    mat mat_data_p=MM_to_mat(filename_p);
    mat_data.append_col(mat_data_p.get_col(0));
  }
  cout<<mat_data<<endl;
  mat_to_MM(mat_data,"data/toric/bp_converge2/toric_S_size_9.mm_result");
}

int bp_converge1(){//result collection for bp decoding //1000 division for p
  //  cout<<"another thing"<<endl;
  //size 5,9,25,35
  char filename_prefix[]="data/toric/bp_converge/toric_S_size_35.mm";
  mat mat_data;
  for (double p=0.001;p<0.050;p+=0.001){
    char filename_p[255];
    sprintf(filename_p,"%s%.3f",filename_prefix,p);
    cout<<filename_p<<endl;
    mat mat_data_p=MM_to_mat(filename_p);
    mat_data.append_col(mat_data_p.get_col(0));
  }
  cout<<mat_data<<endl;
  mat_to_MM(mat_data,"data/toric/bp_converge/toric_S_size_35.mm_result");
}




//collect some individual results and put them into a single file.
int run1(){//for 25 and 9 qubit code
 
  char filename_prefix[]="data/bacon/error_and_failure/BaconShor_size_25_p.mm";//for 25 qubit code
  //  char filename_prefix[]="data/bacon/error_and_failure/BaconShor_size_9_p.mm";//for 9 qubit code
  mat mat_data;  
  for (double p=0.010;p<0.050;p+=0.003){
    char filename_p[255];
    sprintf(filename_p,"%s%.3f",filename_prefix,p);
    mat mat_data_p=MM_to_mat(filename_p);
    mat_data.append_col(mat_data_p.get_col(0));
    
  }
  cout<<mat_data<<endl;
  cout<<"done"<<endl;
  mat_to_MM(mat_data,"data/bacon/error_and_failure/BaconShor_size_25_p.mm");
  return 0;
  
}
int run2(){
  char filename_prefix[]="data/bacon/error_and_failure/BaconShor_size_51_p.mm";
  //run 0-12
  int last_run=19;//start from 1
  mat mat_data;
  for(int i=1;i<last_run+1;i++){
    
  
    mat mat_data_run;  
    for (double p=0.010;p<0.050;p+=0.003){
      char filename_p[255];
      sprintf(filename_p,"%s%d%.3f",filename_prefix,i,p);
      mat mat_data_p=MM_to_mat(filename_p);
      mat_data_run.append_col(mat_data_p.get_col(0));
    
    }
    mat_data.append_row(mat_data_run.get_row(1));
  }
  cout<<mat_data<<endl;
  cout<<"done"<<endl;

  vec average(mat_data.cols());
  average.zeros();
  for(int i=0;i<mat_data.rows();i++){
    average += mat_data.get_row(i);
  }
  average /= 1.0*mat_data.rows();
  cout<<"average"<<endl<<average<<endl;
  mat_data.append_row(average);//the last row gives the average
  mat_to_MM(mat_data,"data/bacon/error_and_failure/BaconShor_size_51_p.mm");
  return 0;
}

int main(){
  //  cout<<"something"<<endl;
  bp_converge2();
  //bp_converge1();
  //run1();
  //run2();
  return 0;
}
