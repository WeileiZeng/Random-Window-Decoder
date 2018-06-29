//Mar 3, 2018 Weilei Zeng
//convert data from .mm file to .gnudat to make it ready for gnuplot
#include "my_lib.h"
#include <itpp/itcomm.h>
#include <stdio.h>
#include <fstream>
using namespace std;
using namespace itpp;


int collect3(char * error_folder, int size){
  //direct rand not included
  //z error only //result collection for bp decoding and rand decoding //100,000 division for p
  //rand directly, BP and rand after BP
  //size 5,9,25,35,..,7,11,13
  
  //  char * error_folder = "data/toric/bp_decoding";
  char filename_prefix[255];
  sprintf(filename_prefix,"%s/toric_S_size_%d.mm_rate",error_folder,size);
  //GF2mat mat_data;
  FILE *fout;//file to save the data
  char filename_out[255];
  sprintf(filename_out,"%s/gnuplot/rate_versus_p_size_%d.gnudat",error_folder,size);
  //cout<<filename_prefix<<endl;
  //  cout<<filename_out<<endl;
  fout =fopen(filename_out,"a");
  
  for (double p=0.00100;p<0.05001;p+=0.00100){

    char filename_p[255];
    sprintf(filename_p,"%s%.5f",filename_prefix,p);
    //    cout<<filename_input<<endl;
    //GF2mat E_input=MM_to_GF2mat(filename_p);
    //    mat_data.append_col(mat_data_p.get_col(0));
    
    //    char filename_input[255];
    //    sprintf(filename_input,"%s_input",filename_p);
    //    GF2mat E_input = MM_to_GF2mat(filename_input);
    GF2mat E_input = get_GF2mat(filename_p,"_input");

    /*
    char filename_input_output_bad[255];
    sprintf(filename_input_output_bad,"%s_input_output_bad",filename_p);
    GF2mat E_input_output_bad = MM_to_GF2mat(filename_input_output_bad);
    */
    /*  
    char filename_output_nonconverge[255];
    sprintf(filename_output_nonconverge,"%s_output_nonconverge",filename_p);
    GF2mat E_output_nonconverge = MM_to_GF2mat(filename_output_nonconverge);
    
    char filename_output_nonconverge_output_bad[255];
    sprintf(filename_output_nonconverge_output_bad,"%s_output_bad",filename_output_nonconverge);
    GF2mat E_output_nonconverge_output_bad = MM_to_GF2mat(filename_output_nonconverge_output_bad);
    */
    GF2mat E_output_nonconverge = get_GF2mat(filename_p,"_output_nonconverge");
    GF2mat E_output_nonconverge_output_bad = get_GF2mat(filename_p,"_output_nonconverge_output_bad");
    
    double N=E_input.rows()-1.0;
    //double N=E_output_converge.rows()+E_output_nonconverge.rows()-2.0;
    //minus one cause the first raw is an extra zero vector.
    //calculate failure rate
    double P_nc = (E_output_nonconverge.rows()-1.0)/N;
    double P_nc_bad = (E_output_nonconverge_output_bad.rows()-1.0)/N;
    //double P_bad = (E_input_output_bad.rows()-1.0)/N;  
    // cout<<"N="<<N<<",P_nc = "<<P_nc<<",P_bad = "<<P_bad<<endl;
    //calculate average weight   
    // double w_i=get_error_density(E_input);   //average weight of input error, appx = p
    double w_temp = get_error_density(E_output_nonconverge);
    /*    if (E_output_nonconverge.rows()==1){
      w_temp=0;
    }else{      
      w_temp = E_output_nonconverge.get_submatrix(1,0,E_output_nonconverge.rows()-1,E_output_nonconverge.cols()-1).density();//average weight of output errors in the non convergent cases
      }*/
    double w_nc = w_temp*P_nc;//average weight of the output errors in BP decoding
    fprintf(fout,"%f\t%f\t%f\t%f\t%f\n",	                                    p, P_nc,P_nc_bad, w_temp,w_nc);  
  }
  fclose(fout);
  cout<<"done converting data for size "<<size<<endl;
  return 0;
}




/*
int collect2new(char * error_folder, int size){//z error only //result collection for bp decoding and rand decoding //100,000 division for p
  //rand directly, BP and rand after BP
  //size 5,9,25,35,..,7,11,13

  //  char * error_folder = "data/toric/bp_decoding";
  char filename_prefix[255];
  sprintf(filename_prefix,"%s/toric_S_size_%d.mm_rate",error_folder,size);
  //GF2mat mat_data;
  FILE *fout;//file to save the data
  char filename_out[255];
  sprintf(filename_out,"%s/gnuplot/rate_versus_p_size_%d.gnudat",error_folder,size);
  fout =fopen(filename_out,"a");
  
  for (double p=0.00100;p<0.05001;p+=0.00100){

    char filename_p[255];
    sprintf(filename_p,"%s%.5f",filename_prefix,p);

    
    double N=E_input.rows()-1.0;    //minus one cause the first raw is an extra zero vector.
    
    //calculate failure rate
    double P_bad = (E_output_bad.rows()-1.0)/N;  
    // cout<<"N="<<N<<",P_nc = "<<P_nc<<",P_bad = "<<P_bad<<endl;

    //calculate average weight   
    double w_bad_temp = get_error_density(E_output_bad);
    doble w_bad = w_bad_temp * P_bad;

    
    fprintf(fout,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
	    p, P_nc,P_nc_bad,P_bad, w_i,w_nc_temp,w_nc, w_nc_bad_temp,w_nc_bad, w_bad_temp,w_bad);    
  }
  fclose(fout);
  cout<<"finish converting data for size "<<size<<endl;
  return 0;
}
*/

int collect_perm(char *filename_output, char * filename_output_bad, char * filename_gnudat, int perm){
  //z error only
  //rate (weight) versus # of permutation
  //  char * error_folder = */p1000;

  FILE *fout;//gnudat file to save the data
  fout =fopen(filename_gnudat,"a");
  //  fprintf(fout,"#rate (weight) versus # of permutation\n	#perm\tP_bad\t\tw_bad_temp\tw_bad\n");
  GF2mat E_output,E_output_bad;

  //  for (int perm=1;perm<21;perm++){
    E_output     = MM_to_GF2mat(filename_output);
    E_output_bad = MM_to_GF2mat(filename_output_bad);
    // E_output_bad=get_GF2mat(error_folder,"perm",perm,filename_output_bad);
    // E_output=get_GF2mat(error_folder,"perm",perm,filename_output);
    double N=E_output.rows()-1.0;    //minus one cause the first raw is an extra zero vector.    
    //calculate failure rate
    double P_bad = (E_output_bad.rows()-1.0)/N;  
    // cout<<"N="<<N<<",P_nc = "<<P_nc<<",P_bad = "<<P_bad<<endl;

    //calculate average weight   
    double w_bad_temp = get_error_density(E_output_bad);
    double w_bad = w_bad_temp * P_bad;

    fprintf(fout,"%d\t%f\t%f\t%f\n",
	    perm*5, P_bad, w_bad_temp,w_bad);    
    // }
  fclose(fout);
  cout<<"get data for perm check when perm = "<<perm<<endl;
  return 0;
}


int collect2(char * error_folder, int size){//z error only //result collection for bp decoding and rand decoding //100,000 division for p
  //rand directly, BP and rand after BP
  //size 5,9,25,35,..,7,11,13

  //  char * error_folder = "data/toric/bp_decoding";
  char filename_prefix[255];
  sprintf(filename_prefix,"%s/toric_S_size_%d.mm_rate",error_folder,size);
  //GF2mat mat_data;
  FILE *fout;//file to save the data
  char filename_out[255];
  sprintf(filename_out,"%s/gnuplot/rate_versus_p_size_%d.gnudat",error_folder,size);
  //cout<<filename_prefix<<endl;
  //  cout<<filename_out<<endl;
  fout =fopen(filename_out,"a");
  
  for (double p=0.00100;p<0.05001;p+=0.00100){

    char filename_p[255];
    sprintf(filename_p,"%s%.5f",filename_prefix,p);
    //    cout<<filename_input<<endl;
    //GF2mat E_input=MM_to_GF2mat(filename_p);
    //    mat_data.append_col(mat_data_p.get_col(0));

    char filename_input[255];
    sprintf(filename_input,"%s_input",filename_p);
    GF2mat E_input = MM_to_GF2mat(filename_input);

    char filename_input_output_bad[255];
    sprintf(filename_input_output_bad,"%s_input_output_bad",filename_p);
    GF2mat E_input_output_bad = MM_to_GF2mat(filename_input_output_bad);
    
    /*char filename_output_converge[255];
    sprintf(filename_output_converge,"%s_output_converge",filename_p);
    GF2mat E_output_converge = MM_to_GF2mat(filename_output_converge);*/
    
    char filename_output_nonconverge[255];
    sprintf(filename_output_nonconverge,"%s_output_nonconverge",filename_p);
    GF2mat E_output_nonconverge = MM_to_GF2mat(filename_output_nonconverge);
    
    char filename_output_nonconverge_output_bad[255];
    sprintf(filename_output_nonconverge_output_bad,"%s_output_bad",filename_output_nonconverge);
    GF2mat E_output_nonconverge_output_bad = MM_to_GF2mat(filename_output_nonconverge_output_bad);
    
    
    double N=E_input.rows()-1.0;
    //double N=E_output_converge.rows()+E_output_nonconverge.rows()-2.0;
    //minus one cause the first raw is an extra zero vector.
    //calculate failure rate
    double P_nc = (E_output_nonconverge.rows()-1.0)/N;
    double P_nc_bad = (E_output_nonconverge_output_bad.rows()-1.0)/N;
    double P_bad = (E_input_output_bad.rows()-1.0)/N;  
    // cout<<"N="<<N<<",P_nc = "<<P_nc<<",P_bad = "<<P_bad<<endl;

    //calculate average weight   
    double w_i=get_error_density(E_input);    //average weight of input error, appx = p

    double w_temp=get_error_density(E_output_nonconverge);
    //average weight of output errors in the non convergent cases    
    double w_nc = w_temp*P_nc;//average weight of the output errors in BP decoding
    
    fprintf(fout,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
	    p, P_nc,P_nc_bad,P_bad, w_i,w_temp,w_nc);    
  }
  fclose(fout);
  cout<<"done..."<<endl;
  return 0;
}


int collect1(int size){//result collection for bp decoding and rand decoding //100,000 division for p
  //rand directly, BP and rand after BP
  //size 5,9,25,35,..,7,11,13
  //  int size=7;
  char filename_prefix[255];
  sprintf(filename_prefix,"data/toric/bp_converge3/toric_S_size_%d.mm_rate",size);
  //GF2mat mat_data;
  FILE *fout;//file to save the data
  char filename_out[255];
  sprintf(filename_out,"data/toric/bp_converge3/gnuplot/rate_versus_p_size_%d.gnudat",size);
  //cout<<filename_prefix<<endl;
  //  cout<<filename_out<<endl;
  fout =fopen(filename_out,"a");
  
  for (double p=0.00100;p<0.05001;p+=0.00100){

    char filename_p[255];
    sprintf(filename_p,"%s%.5f",filename_prefix,p);
    //    cout<<filename_input<<endl;
    //GF2mat E_input=MM_to_GF2mat(filename_p);
    //    mat_data.append_col(mat_data_p.get_col(0));

    char filename_input[255];
    sprintf(filename_input,"%s_input",filename_p);
    GF2mat E_input = MM_to_GF2mat(filename_input);

    char filename_input_output_bad[255];
    sprintf(filename_input_output_bad,"%s_input_output_bad",filename_p);
    GF2mat E_input_output_bad = MM_to_GF2mat(filename_input_output_bad);
    
    /*char filename_output_converge[255];
    sprintf(filename_output_converge,"%s_output_converge",filename_p);
    GF2mat E_output_converge = MM_to_GF2mat(filename_output_converge);*/
    
    char filename_output_nonconverge[255];
    sprintf(filename_output_nonconverge,"%s_output_nonconverge",filename_p);
    GF2mat E_output_nonconverge = MM_to_GF2mat(filename_output_nonconverge);
    
    char filename_output_nonconverge_output_bad[255];
    sprintf(filename_output_nonconverge_output_bad,"%s_output_bad",filename_output_nonconverge);
    GF2mat E_output_nonconverge_output_bad = MM_to_GF2mat(filename_output_nonconverge_output_bad);
    
    
    double N=E_input.rows()-1.0;
    //double N=E_output_converge.rows()+E_output_nonconverge.rows()-2.0;
    //minus one cause the first raw is an extra zero vector.
    //calculate failure rate
    double P_nc = (E_output_nonconverge.rows()-1.0)/N;
    double P_nc_bad = (E_output_nonconverge_output_bad.rows()-1.0)/N;
    double P_bad = (E_input_output_bad.rows()-1.0)/N;  
    // cout<<"N="<<N<<",P_nc = "<<P_nc<<",P_bad = "<<P_bad<<endl;

    //calculate average weight   
    double w_i=E_input.get_submatrix(1,0,E_input.rows()-1,E_input.cols()-1).density();//average weight of input error, appx = p
    double w_temp = E_output_nonconverge.get_submatrix(1,0,E_output_nonconverge.rows()-1,E_output_nonconverge.cols()-1).density();//average weight of output errors in the non convergent cases
    double w_nc = w_temp*P_nc;//average weight of the output errors in BP decoding

    //break;
    fprintf(fout,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
	    p, P_nc,P_nc_bad,P_bad, w_i,w_temp,w_nc);
    
  }
  fclose(fout);
  //cout<<mat_data<<endl;
  //mat_to_MM(mat_data,"data/toric/bp_converge2/toric_S_size_9.mm_result");
  cout<<"done..."<<endl;
  return 0;
}






int bp_converge3(int size){//result collection for bp decoding and rand decoding //100,000 division for p
  //cout<<"another thing"<<endl;
  //size 5,9,25,35,..,7,11,13
  //  int size=7;
  char filename_prefix[255];
  sprintf(filename_prefix,"data/toric/bp_converge3/toric_S_size_%d.mm_rate",size);
  //GF2mat mat_data;
  FILE *fout;//file to save the data
  char filename_out[255];
  sprintf(filename_out,"data/toric/bp_converge3/gnuplot/rate_versus_p_size_%d.gnudat",size);
  //cout<<filename_prefix<<endl;
  //  cout<<filename_out<<endl;
  fout =fopen(filename_out,"a");
  
  for (double p=0.00100;p<0.05001;p+=0.00100){

    char filename_p[255];
    sprintf(filename_p,"%s%.5f",filename_prefix,p);
    //    cout<<filename_input<<endl;
    //GF2mat E_input=MM_to_GF2mat(filename_p);
    //    mat_data.append_col(mat_data_p.get_col(0));
    char filename_output_converge[255];
    sprintf(filename_output_converge,"%s_output_converge",filename_p);
    GF2mat E_output_converge = MM_to_GF2mat(filename_output_converge);
    char filename_output_nonconverge[255];
    sprintf(filename_output_nonconverge,"%s_output_nonconverge",filename_p);
    GF2mat E_output_nonconverge = MM_to_GF2mat(filename_output_nonconverge);
    
    char filename_output_nonconverge_output_bad[255];
    sprintf(filename_output_nonconverge_output_bad,"%s_output_bad",filename_output_nonconverge);
    GF2mat E_output_nonconverge_output_bad = MM_to_GF2mat(filename_output_nonconverge_output_bad);
    
    
    double N=E_output_converge.rows()+E_output_nonconverge.rows()-2.0;
    //minus one cause the first raw is an extra zero vector.
    double P_nc = (E_output_nonconverge.rows()-1.0)/N;
    double P_bad = (E_output_nonconverge_output_bad.rows()-1.0)/N;
    //cout<<"N="<<N<<",P_nc = "<<P_nc<<",P_bad = "<<P_bad<<endl;
    fprintf(fout,"%f\t%f\t%f\n",p,P_nc,P_bad);
    
  }
  fclose(fout);
  //cout<<mat_data<<endl;
  //mat_to_MM(mat_data,"data/toric/bp_converge2/toric_S_size_9.mm_result");
  cout<<"done..."<<endl;
  return 0;
}



int merge_file(char * error_folder,int size){
  //merge two files into one file;
  //E_input_converge and E_input_nonconverge into E_input

  //size 5,9,25,35,..,7,11,13
  //  int size=7;
  //  char * error_folder = "data/toric/bp_decoding";
  char filename_prefix[255];
  sprintf(filename_prefix,"%s/toric_S_size_%d.mm_rate",error_folder,size);
  
  for (double p=0.00100;p<0.05001;p+=0.00100){//use 0.05001 cause 0.05000 is included

    char filename_p[255];
    sprintf(filename_p,"%s%.5f",filename_prefix,p);
    char filename_input_converge[255];
    sprintf(filename_input_converge,"%s_input_converge",filename_p);
    GF2mat E_input_converge = MM_to_GF2mat(filename_input_converge);
    char filename_input_nonconverge[255];
    sprintf(filename_input_nonconverge,"%s_input_nonconverge",filename_p);
    GF2mat E_input_nonconverge = MM_to_GF2mat(filename_input_nonconverge);
    GF2mat E_input;
    if (E_input_nonconverge.rows()==1){
      E_input=E_input_converge;
    }else if (E_input_converge.rows()==1){
      E_input=E_input_nonconverge;
    }else{
      E_input = E_input_converge.concatenate_vertical(				      E_input_nonconverge.get_submatrix(1,0,E_input_nonconverge.rows()-1,E_input_nonconverge.cols()-1)	        );
    }
   
    char filename_input[255];
    sprintf(filename_input,"%s_input",filename_p);
    GF2mat_to_MM(E_input,filename_input);
  }
  cout<<"finish merge files for size "<<size<<endl;
  return 0;
}

int main(int argc, char ** argv){
  if (argc<5 ){ 
    //cout<<"Please enter the error_folder and size of the code..."<<endl;
    cout<<"Goodbye"<<endl;
  }
  //  char * error_folder = argv[1];
  //  int size =atof(argv[2]);
  //  cout<<size<<endl;
  //bp_converge3(size);
  //    merge_file(error_folder,size);//merge two files into one
  //  collect2(error_folder,size);//BP, rand after bp and rand
  //collect2new(error_folder,size);//BP, rand after bp and rand
  // collect3(error_folder,size);//rand directly not included
  //collect1(size);

  //  char * error_folder = argv[1];
  char *filename_output = argv[1];
  char * filename_output_bad = argv[2];
  char * filename_gnudat = argv[3];
  int perm = atof(argv[4]);
  collect_perm( filename_output, filename_output_bad,  filename_gnudat, perm);

  return 0;
}
