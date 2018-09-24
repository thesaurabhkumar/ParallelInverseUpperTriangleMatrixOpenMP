#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void printMatrixKonsole(double** Mat, int x, int y)
{
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            printf("%4.2f ", Mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");
}

void printMatrix(double** Mat, int x, int y, FILE *fp)
{
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            //printf("%4.2f ", Mat[i][j]);
            fprintf(fp, "%4.5lf ", Mat[i][j]);
        }
        //printf("\n");
        fprintf(fp, "\n");
    }
    //printf("\n");
    //printf("\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");
}

double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}

/*
   Find the cofactor matrix of a square matrix
*/
void CoFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}

/*
   Transpose of a square matrix, do it in place
*/
void Transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}

void initMatrix(double** Mat, int x, int y)
{
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            Mat[i][j] = rand()%2000 + 1;
        }
    }
}

double** upperTrian(double** Mat, int n){
    int k=0;
    for(int i=0; i<n; i++){
        for(int j=0; j<k; j++){
            Mat[i][j] = 0;
            //printf("saurabh\n");
        }
        k++;
    }
    return(Mat);
}

double norm(double** MatA, double** MatB, int n){
    double norm = 0;
    double** mul = matrixMul(MatA, MatB, n);
    //printMatrixKonsole(MatA, n, n);
    //printMatrixKonsole(MatB, n, n);
    //printf("mul\n");
    //printMatrixKonsole(mul, n, n);
    //printf("n=%d\n", n);
    for(int i=0; i<n; i++){
	for(int j=0; j<n; j++){
	    if(i==j){
                mul[i][j] = mul[i][j] - 1;
		//printf("sayravg\n");
	    }
	}
    }
    //printf("mul\n");
    //printMatrixKonsole(mul, n, n);

    for(int i=0; i<n; i++){
	for(int j=0; j<n; j++){
	    double prev_norm = norm;
    	    norm = norm + (mul[i][j] * mul[i][j]);
	    if((norm - prev_norm) > 100){
	    	printf("prev_norm=%lf norm=%lf\n", prev_norm, norm);
	    //if(mul[i][j] > 0.0){
	    	printf("i=%d j=%d mul=%lf\n", i,j,mul[i][j]);
	    }
	    //printf("%4.2d  ", mul[i][j]);
	}
	//printf("\n");
    }
    return(norm);
}

int main(int argc, char** argv){
    FILE *fp = fopen("output.txt", "w");
    int n = 0;
    if(argc != 2){
        printf("Incorrect number of arguments: %d. Please Check!\n", argc);
        exit(-1);
    }
    n = atoi(argv[1]);
    double** Mat = allocMatrix(n, n);
    initMatrix(Mat, n, n);
    Mat = upperTrian(Mat,n);
    printMatrix(Mat, n, n, fp);
    
    //double start=omp_get_wtime();
    //printf("compute_inverse time = %10.4e\n", omp_get_wtime() - start);
    //printMatrix(Mat, n, n, fp);
    printMatrix(Inv, n, n, fp);

    printf("\ncomputation norm = %lf\n", norm(Mat, Inv, n));
    fclose(fp);
    return(0);
}


