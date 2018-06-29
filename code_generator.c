#include "mm_write.h"
#include <itpp/itbase.h>
#include <stdio.h>
using namespace itpp;
using namespace std;

GF2mat get_tilde(GF2mat G){
  GF2mat G_t=G.get_submatrix(0,G.cols()/2,G.rows()-1,G.cols()-1).concatenate_horizontal( G.get_submatrix(0,0,G.rows()-1,G.cols()/2-1) );
  return G_t;
}

int bacon_shor(int n, char* filename);

int bravyi(int n, double density,
	   char* filename_A="bravyi_A.dat",char* filename_G="bravyi_G.dat");
//size n x n; density between 0 and 1

//int toric(int L);
int toric(int L, char * filename_S);
int toric_S_x(int L, char * filename_S);
int toric_S_z(int L, char * filename_S);

int generate();

int bs_generate(){
  //readme.txt
  bacon_shor(9,"data/bacon/BaconShor_size_9_G.dat");
  bacon_shor(25,"data/bacon/BaconShor_size_25_G.dat");
  bacon_shor(51,"data/bacon/BaconShor_size_51_G.dat");
  bacon_shor(125,"data/bacon/BaconShor_size_125_G.dat");
  return 0;
}

int main(int args, char ** argv){
  //  Real_Timer t;
  //  t.tic();
  //bs_generate();
  //  bravyi(5,1.0);
  //  generate();
  //bacon_shor(27);
  // toric(75);
  //toric_S_x( atof(argv[1]),argv[2]  );
  toric_S_z( atof(argv[1]),argv[2]  );
  //    toric( atof(argv[1]),argv[2]  );
  //  toric(5);
  // t.toc_print();
  return 0;
}

int generate(){
  //generate some bravyi code
  //step 1: size n 5,7,9,...,21;finished
  //density: 0.6,0.7,0.8,0.9
  //trials 10 for each parameter set
  //Note:  number of zeros in A, n0=n*n*(1-density), number of matrices: 2^(n0+1)
  //step2: size n=(23,33,2) finished
  //step3: size n=(33,43,2)
  for (int n=33;n<43;n+=2){//size
    for (double d=0.5;d<0.91;d+=0.1){//density
      for (int i=0;i<10;i++){//10 for each
	char filename_A[66],filename_G[66];
	sprintf(filename_A,"data/bravyi_%s_size_%d_density_%.2f_v_%d.dat","A",n,d,i);
	sprintf(filename_G,"data/bravyi_%s_size_%d_density_%.2f_v_%d.dat","G",n,d,i);
	//	sprintf(filename_G,"%s%s%s%d%s%.2fv%d%s","data/bravyi_","G","_size_",n,"_density_",d,i,".dat");
	bravyi(n,d,filename_A,filename_G);
      }
    }
  }
  return 0;
}

//can combine S_x and S_z by adding ann extra variable
int toric_S_z(int L, char * filename_S){
  //only have X vertex operators and only check z errors
  //size n=L x L. total spins 2n. toric code is a stabilzier code
  // cout<<"hello here"<<endl;
  GF2mat G((L*L-1),L*L*2);//G is stabilizer matrix for toric code. //all zero when initialized
  //positions of spins/bonds: (horizontal links for X, vertical links for X)
  int row=-1;//the current number of stbailzier generators has been added into G
  //  cout<<"debug :"<<endl;
  for (int i =0;i<L;i++){//(i,j) is the physical position of the qubit. Not the position in G
        // cout<<"\t i="<<i;
    for (int j=0;j<L;j++){
          int pos=i*L+j;
      row++;//row is the number of generators, it is the row in generating matrix
      /*
      //vertex for X
      //horizontal links
      G.set(row,(pos+L)%(L*L),1);
      G.set(row,(pos+L+1 -(j+1)/L*L )%(L*L),1);
      //vertical links
      G.set(row,L*L+pos+1 -(j+1)/L*L,1);
      G.set(row,L*L+(pos+L+1 -(j+1)/L*L )%(L*L),1);
      */

      //plaquette for Z
      //horizontal links
      G.set(row,+pos,1);
      G.set(row,+(pos+L) % (L*L)   ,1);
      //vertical links
      G.set(row,L*L+pos,1);
      G.set(row,L*L+ pos+1 -  (j+1)/L*L  ,1);// -(j+1)/L*L  for periodic boundary condition
      
      if (j==L-2 && i==L-1){
	break;//remove the last one.
      }
    }
  }
  cout<<"finish generating S_z (plaquete with Z operators) for lattice size "<<L<<" x "<<L<<endl;
  // cout<<"Size of G is "<<G.rows()<<" x "<<G.cols()<<endl;
   GF2mat_to_MM(G,filename_S);
  //check
  //  cout<<"G.row_rank() = "<<G.row_rank()<<endl;
  //cout<<G<<endl;
  return 0;
}

int toric_S_x(int L, char * filename_S){
  //only have X vertex operators and only check z errors
  //size n=L x L. total spins 2n. toric code is a stabilzier code
  GF2mat G((L*L-1),L*L*2);//G is stabilizer matrix for toric code. //all zero when initialized
  //positions of spins/bonds: (horizontal links for X, vertical links for X)
  int row=-1;//the current number of stbailzier generators has been added into G
  for (int i =0;i<L;i++){//(i,j) is the physical position of the qubit. Not the position in G
    //    cout<<"\t i="<<i;
    for (int j=0;j<L;j++){
      int pos=i*L+j;
      row++;//row is the number of generators, it is the row in generating matrix
      //vertex for X
      //horizontal links
      G.set(row,(pos+L)%(L*L),1);
      G.set(row,(pos+L+1 -(j+1)/L*L )%(L*L),1);
      //vertical links
      G.set(row,L*L+pos+1 -(j+1)/L*L,1);
      G.set(row,L*L+(pos+L+1 -(j+1)/L*L )%(L*L),1);
      /*
      //plaquette for Z
      //horizontal links
      G.set(row+L*L-1,2*L*L+pos,1);
      G.set(row+L*L-1,2*L*L+(pos+L) % (L*L)   ,1);
      //vertical links
      G.set(row+L*L-1,3*L*L+pos,1);
      G.set(row+L*L-1,3*L*L+ pos+1 -  (j+1)/L*L  ,1);// -(j+1)/L*L  for periodic boundary condition
      */
      if (j==L-2 && i==L-1){
	break;//remove the last one.
      }
    }
  }
  cout<<"finish generating G for lattice size "<<L<<" x "<<L<<endl;
  // cout<<"Size of G is "<<G.rows()<<" x "<<G.cols()<<endl;
   GF2mat_to_MM(G,filename_S);
  //check
  //  cout<<"G.row_rank() = "<<G.row_rank()<<endl;
  //cout<<G<<endl;
  return 0;
}


int toric(int L, char * filename_S){
  //it takes 5 secs to construct the code with size 200 x 200
  //size n=L x L. total spins 2n. toric code is a stabilzier code
  GF2mat G(2*(L*L-1),L*L*4);//G is stabilizer matrix for toric code. //all zero when initialized
  //positions of spins: (horizontal links for X, vertical links for X,horizontal links for Z, vertical links for Z, )
  int row=-1;//the current number of stbailzier generators hass been added into G
  for (int i =0;i<L;i++){//(i,j) is the physical position of the qubit. Not the position in G
    cout<<"\t i="<<i;
    for (int j=0;j<L;j++){
      int pos=i*L+j;
      row++;//row is the number of generators, it is the row in generating matrix
      //vertex for X
      //horizontal links
      G.set(row,(pos+L)%(L*L),1);
      G.set(row,(pos+L+1 -(j+1)/L*L )%(L*L),1);
      //vertical links
      G.set(row,L*L+pos+1 -(j+1)/L*L,1);
      G.set(row,L*L+(pos+L+1 -(j+1)/L*L )%(L*L),1);

      //plaquette for Z
      //horizontal links
      G.set(row+L*L-1,2*L*L+pos,1);
      G.set(row+L*L-1,2*L*L+(pos+L) % (L*L)   ,1);
      //vertical links
      G.set(row+L*L-1,3*L*L+pos,1);
      G.set(row+L*L-1,3*L*L+ pos+1 -  (j+1)/L*L  ,1);// -(j+1)/L*L  for periodic boundary condition
  
      if (j==L-2 && i==L-1){
	break;//remove the last one.
      }
    }
  }
  cout<<"finish generating G"<<endl;
  GF2mat_to_MM(G,filename_S);
  
  /* GF2mat G_t=get_tilde(G);
  //  cout<<G_t<<endl;
  int weight=0;
  cout<<"check commutation"<<endl;
  GF2mat F=G_t*G.transpose();//this procedure takes too long;
  t.toc_print();
  for (int i =0;i<F.rows();i++){
    cout<<"i= "<<i<<endl;
    for ( int j=0;j< F.cols();j++){
      weight =  F.get(i,j) ?  weight +1 : weight;
    }
  }
  cout<<"weight of F is "<<weight<<endl;
  */  
  
  return 0;
}



int bacon_shor(int n,char* filename){//size n x n
  int s=n*n;//size s
  GF2mat G( 2*n*(n-1),s*2 );
  int x,y;
  for (int i=0;i<n;i++){
    for(int j=0;j<n-1;j++){
      //pairs of X operators in the same row
      x=i*(n-1)+j;
      y=x+i;
      G.set(x,y,1);
      G.set(x, y+1,1);
      //pairs of Z operators in the same column
      x=n*(n-1)+j+i*(n-1);
      y=s+j*n+i;
      G.set(x,y,1);
      G.set(x,y+n,1);

    }
  }
  //    cout<<G<<endl;
    GF2mat_to_MM(G,filename);//"G_BS_size11.dat");
  return 0;
}


int bravyi(int n, double density,
	   char* filename_A,	   char* filename_G){//size n x n; density between 0 and 1
  RNG_randomize();
  GF2mat A(n,n);
  for (int i=0;i<n;i++){
    for ( int j=0;j<n;j++){
      if(randu()<density){
	A.set(i,j,1);
      }
    }
  }
  //  cout<<"A="<<A<<endl;
  GF2mat_to_MM(A,filename_A);
  //construct G from A
  GF2mat G( 2*n*(n-1),2*n*n );//big redundant zero matrix. fully used in the case of Bacon Shor code
  int num=0;//number of generator added to G
  //X generators, pairs of X in same row
  for (int i=0;i<n;i++){
    int a=A.get(i,0);
    int single=a;//if the row has exact one element =1 .
    int last=a ;//if there exists an 1 in previous element of this row
    int last_pos=0;
    for (int j=1;j<n;j++){
      if(A.get(i,j) ){//act for 1 and pass for 0
	if(last){
	  //create a pair of generator
	  G.set(num,i*n+last_pos,1);
	  G.set(num,i*n+j,1);
	  num++;
	  single=false;
	}else{
	  //record that we found an 1
	  last=true;
	  single=true;
	}
	last_pos=j;//record this position for future use
      }
    }
    if(single){
      //create a single operator
      G.set(num,i*n+last_pos,1);
      num++;
    }
  }

  //Z generators, pairs of X in same column
  for (int i=0;i<n;i++){//i for columns
    int a=A.get(0,i);
    int single=a;//if the column has exact one element =1 .
    int last=a ;//if there exists an 1 in previous element of this column
    int last_pos=0;
    for (int j=1;j<n;j++){//j for rows
      if(A.get(j,i) ){//act for 1 and pass for 0
	if(last){
	  //create a pair of generator
	  G.set(num,n*n+last_pos*n+i,1);
	  G.set(num,n*n+j*n+i,1);
	  num++;
	  single=false;
	}else{
	  //record that we found an 1
	  last=true;
	  single=true;
	}
	last_pos=j;//record this position for future use
      }
    }
    if(single){
      //create a single operator
      G.set(num,n*n+last_pos*n+i,1);
      num++;
    }
  }
  G.set_size(num,G.cols(),true);
  //cout<<"G="<<G<<endl;
  //remove zero columns in G
  int zero=0;//
  for (int i=0;i<n;i++){
    for ( int j=0;j<n;j++){
      if (A.get(i,j) ){//1: keep the column and move it to the right place.
	G.swap_cols(i*n+j,i*n+j-zero);
	G.swap_cols(n*n+i*n+j,n*n+i*n+j-zero);
      }else{//0:  pass
	zero++;
	
      }
	
    }
  }
  int wt=n*n-zero;//weight of A
  G=(  G.get_submatrix(0,0,G.rows()-1,wt-1)
       ).concatenate_horizontal(
				G.get_submatrix(  0,n*n,G.rows()-1,n*n+wt-1  )   );

  //cout<<"G="<<G<<endl;
  GF2mat_to_MM(G,filename_G);
  
  return 0;
}
