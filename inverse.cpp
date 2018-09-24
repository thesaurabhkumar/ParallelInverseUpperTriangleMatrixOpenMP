#include <iostream>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>


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

double** allocMatrix(int x, int y){
    double **matrix = NULL;
    matrix = (double**)malloc(x * sizeof(double*));
    for(int i=0; i<x; i++){
        matrix[i] = (double*)malloc(y * sizeof(double));
    }
    return(matrix);
}

double** matrixMul(double** Mat1, double** Mat2, int n)
{
    double** res = allocMatrix(n, n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            //res[i][j] = 0;
            for(int k=0; k<n; k++){
                //printf("res=%lf  %lf  %lf\n", res[i][j], Mat1[i][k], Mat2[k][j]);
                res[i][j] += Mat1[i][k]*Mat2[k][j];
                //printf("i=%d j=%d k=%d res=%lf\n", i, j, k, res[i][j]);
            }
        }
    }
    //printf("A\n");
    //printMatrixKonsole(Mat1, n, n);
    //printf("B\n");
    //printMatrixKonsole(Mat2, n, n);
    //printf("Res\n");
    //printMatrixKonsole(res, n, n);
    return(res);
}

double** seqInv(double **Mat, int x1, int y1, int x2, int y2){
    double **Inv = allocMatrix(x2 - x1 + 1, y2 - y1 + 1);
    double det = Mat[x1][y1]*Mat[x2][y2] - Mat[x1][y2]*Mat[x2][y1];
    Inv[0][0] = Mat[x2][y2]/det;
    Inv[0][1] = -Mat[x1][y2]/det;
    Inv[1][0] = -Mat[x2][y1]/det;
    Inv[1][1] = Mat[x1][y1]/det;
    //printf("x1=%d y1=%d x2=%d y2=%d\n", x1, y1, x2, y2);
	//printMatrix(Inv, x2-x1+1, y2-y1+1);
    return(Inv);
}

double** sliceMatrix(double** Mat, int x1, int y1, int x2, int y2)
{
    double** slice = allocMatrix(x2-x1+1, y2-y1+1);
    for(int i=0; i<x2-x1+1; i++){
        for(int j=0; j<y2-y1+1; j++){
            slice[i][j] = Mat[x1+i][y1+j];
        }
    }
    //printf("SLICE x1=%d y1=%d x2=%d y2=%d\n", x1, y1, x2, y2);
	//printMatrix(slice, x2-x1+1, y2-y1+1);
    return(slice);
}

double** negMatrix(double** Mat, int x, int y){
    double** res = allocMatrix(x+1, y+1);
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            res[i][j] = -Mat[i][j];
        }
    }
    return(res);
}

void compute_inverse(double** Mat, int x1, int y1, int x2, int y2)
{
    int x = x2 - x1 + 1;
    int y = y2 - y1 + 1;
    int n = x;

    if(x != y){
        printf("Trying to calculate inverse of a matrix of order %d by %d\n", x, y);
        exit(-1);
    }
    if(n==1) {
        Mat[x1][y1] = 1/Mat[x1][y1];
        return ;
    }

    int n1=(n-1)/2;

	#pragma omp parallel if(n>=200)
    {
        #pragma omp single
        {
            #pragma omp task default(none) shared(Mat,x1,y1,n1)
            compute_inverse(Mat,x1,y1,x1+n1,y1+n1);
            #pragma omp task default(none) shared(Mat,x1,y1,x2,y2,n1)
            compute_inverse(Mat,x1+n1+1,y1+n1+1,x2,y2);
            #pragma omp taskwait
        }
    }

    #pragma omp parallel default(none) shared(Mat,x1,y1,x2,y1,n)
    {
        double *mRes = (double*)malloc((n1+1)*sizeof(double));
        #pragma omp for
        for(int j=y1+n1+1; j<y2+1; ++j){
            for(int i=x1; i<x1+n1+1; ++i){
                mRes[i-x1] = 0;
                for(int k=0; k<n1+1; ++k){
                    mRes[i-x1] += -1*Mat[i][y1+k]*Mat[x1+k][j];
            }
        }
        for(int l=0;l<n1+1;++l){
            Mat[x1+l][j]=mRes[l];
        }
    }

    double *mRes2 = (double*)malloc((n-n1-1)*sizeof(double));

    #pragma omp for
        for(int i=x1;i<x1+n1+1;++i){
            for(int j=y1+n1+1;j<y1+n;++j){
                mRes2[j-y1-n1-1] = 0;
                for(int k=0;k<n-n1-1;++k){
                    mRes2[j-y1-n1-1] += Mat[i][y1+n1+1+k]*Mat[x1+n1+1+k][j];
                }
            }
            for(int l=0;l<n-n1-1;++l){
                Mat[i][y1+n1+1+l] = mRes2[l];
            }

        }
    }

}

double** upperTrian(double** Mat, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<i; j++){
            Mat[i][j] = 0;
        }
    }
    return(Mat);
}

void initMatrix(double** Mat, int x, int y)
{
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            Mat[i][j] = rand()%2000 + 1;
        }
    }
    for(int j=0; j<y; j++){
        int sum = 0;
        for(int i=0; i<x; i++){
            sum += Mat[i][j];
        }
        Mat[j][j] += sum;
    }
}

void copyMatrix(double** inMat, double** outMat, int n)
{
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            outMat[i][j] = inMat[i][j];
        }
    }

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
    	    norm = norm + (mul[i][j] * mul[i][j]);
	    }
    }
    norm = sqrt(norm);
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
    double** inMat = allocMatrix(n, n);
    double** outMat = allocMatrix(n, n);
    initMatrix(inMat, n, n);
    inMat = upperTrian(inMat,n);
    copyMatrix(inMat, outMat, n);
    printMatrix(inMat, n, n, fp);
    
    //double start_time = omp_get_wtime();
    compute_inverse(outMat, 0, 0, n-1, n-1);
    //printf("compute_inverse time = %10.5e\n", omp_get_wtime() - start_time);
    printMatrix(outMat, n, n, fp);

    //start_time = omp_get_wtime();
    printf("\nComputation Norm = %lf\n", norm(inMat, outMat, n));
    //printf("norm calculation time = %10.5e\n", omp_get_wtime() - start_time);
    fclose(fp);
    return(0);
}


