#include <stdio.h>
//#include <bits/stdc++.h>
//#include <omp.h>
#include <vector>

#define v vector
#define REP(i,a,n) for(int i=a;i<n;++i)
#define ld long double

int num_of_threads;
int threads_error_calc;

using namespace std;

void compute_inverse(v<v<double>> &R, int r1, int c1, int r2, int c2){
  
  int n=r2-r1+1;
  
  if(n==1) {    //base case
    
    R[r1][c1] = 1/R[r1][c1];
     
    return ;
  }    //base case ends

  int n1=(n-1)/2;  //index
  
    #pragma omp parallel if(n>=500)
    {
      #pragma omp single
      {
          #pragma omp task default(none) shared(R,r1,c1,n1)
          compute_inverse(R,r1,c1,r1+n1,c1+n1);
          
          #pragma omp task default(none) shared(R,r1,c1,r2,c2,n1)
          compute_inverse(R,r1+n1+1,c1+n1+1,r2,c2);
          
          #pragma omp taskwait
      }
    }

  
    #pragma omp parallel default(none) shared(R,r1,c1,c2,n1,n)
    {
    double mRes[n1+1];

  
    #pragma omp for
    for(int j=c1+n1+1;j<c2+1;++j){
    //REP(j,c1+n1+1,c2+1){
        //REP(i,r1,r1+n1+1){
        for(int i=r1;i<r1+n1+1;++i){
        mRes[i-r1]=0;    //uninitialization caused random errors

            for(int k=0;k<n1+1;++k){
            //REP(k,0,n1+1){
            mRes[i-r1]+=-1*R[i][c1+k]*R[r1+k][j];
            }
        }

        for(int l=0;l<n1+1;++l){
        //REP(l,0,n1+1){
            R[r1+l][j]=mRes[l];
        }
    
    }
  
    double mRes2[n-n1-1];
  
    #pragma omp for
    for(int i=r1;i<r1+n1+1;++i){
    //REP(i,r1,r1+n1+1){
        for(int j=c1+n1+1;j<c1+n;++j){
        //REP(j,c1+n1+1,c1+n){
            mRes2[j-c1-n1-1]=0;
            for(int k=0;k<n-n1-1;++k){
            //REP(k,0,n-n1-1){
                mRes2[j-c1-n1-1]+=R[i][c1+n1+1+k]*R[r1+n1+1+k][j];
            }
        }

        for(int l=0;l<n1-n1-1;++l){
        //REP(l,0,n-n1-1){
            R[i][c1+n1+1+l]=mRes2[l];
        }
    
    }
    }
  
  return;
}


void Rinverse(int &n){
  default_random_engine generator (0);
  uniform_real_distribution<float> distribution (0,1);
  
  v<v<double>> A(n,v<double> (n,0));

  REP(i,0,n) REP(j,0,n) A[i][j]=1-distribution(generator);
  
  double sum=0.0;
  REP(j,0,n) {
    sum=0.0;
        REP(i,0,n)
            sum=sum+A[i][j];
        A[j][j] += sum;
  }
  
  //make upper triangular
  REP(i,0,n) REP(j,0,i) A[i][j]=0.0;
  
  //compute inverse
  v<v<double>> Ri=A ;
  
  double start=omp_get_wtime();
  compute_inverse(Ri,0,0,n-1,n-1);
  printf("Time - compute_inverse:    %10.4e\n", omp_get_wtime()-start);
  
  start=omp_get_wtime();
  //error calc
  long double e=0;
  
  #pragma omp parallel for collapse(2) shared(A,Ri,n) reduction(+:e) num_threads(threads_error_calc)
  REP(i,0,n)
    REP(j,0,n){
      long double t=0;
      REP(k,0,n){
        t=t+A[i][k]*Ri[k][j];
        
      }

      if(i!=j) e+=pow(t,2);
      else e+=pow(t-1,2);

    }

  printf("Time - error calculation:    %10.4e\n", omp_get_wtime()-start);
  //REP(i,0,n) {REP(j,0,n) cout<<setprecision(4)<<A[i][j]<<" "; cout<<endl;} cout<<endl;
  //REP(i,0,n) {REP(j,0,n) cout<<setprecision(4)<<Ri[i][j]<<" "; cout<<endl;} cout<<endl;
  
  cout<<"Error in computing inverse: "<<setprecision(6)<<sqrt(e)<<endl;
  cout<<endl;
  return; 
}

int main(int argc, char **argv){
  
  int n = atoi(argv[1]);
  num_of_threads = atoi(argv[2]);
  threads_error_calc = atoi(argv[3]);
  
  omp_set_num_threads(num_of_threads);
  
  Rinverse(n);
  
  return 0;
}

