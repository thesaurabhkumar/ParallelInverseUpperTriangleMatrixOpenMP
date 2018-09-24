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

double** compute_inverse(double** Mat, int x1, int y1, int x2, int y2)
{
    int x = x2 - x1 + 1;
    int y = y2 - y1 + 1;
    if(x != y){
        printf("Trying to calculate inverse of a matrix of order %d by %d\n", x, y);
        exit(-1);
    }
    if((x < 2) || (y < 2)){
        printf("Trying to calculate inverse of a small matrix of order %d by %d\n", x, y);
        exit(-1);
    }
    int n = x*y;

    if(n == 4){
	    return(seqInv(Mat, x1, y1, x2, y2));
    } else{
        double **Inv1 = NULL;
        double **Inv2 = NULL;
        double **Inv = allocMatrix(x2 - x1 + 1, y2 - y1 + 1);
        double **Inv3 = NULL;
        int xmid = round((x2 + x1 + 1)/2);
        int ymid = round((y2 + y1 + 1)/2);
	
	//#pragma omp parallel if(n>=100)
	{
	//#pragma omp single
	{
	    //#pragma omp task default(none) shared(Mat, Inv1, x1, y1, xmid, ymid)
	Inv1 = compute_inverse(Mat, x1, y1, xmid-1, ymid-1);
	    //#pragma omp task default(none) shared(Mat, Inv2, x2, y2, xmid, ymid)
	Inv2 = compute_inverse(Mat, xmid, ymid, x2, y2);
	//#pragma omp taskwait
	}
	}
        double** slice = sliceMatrix(Mat, x1, ymid, xmid-1, y2);
        //printf("RECUR x1=%d y1=%d xmid=%d ymid=%d x2=%d y2=%d\n", x1, y1, xmid, ymid, x2, y2);
        //printMatrix(Inv1, xmid-x1, ymid-y1);
        //printMatrix(Inv2, x2-xmid+1, y2-ymid+1);

        //printf("saurabh\n");
	    double **temp1 = matrixMul(Inv1, slice, x2-xmid+1);
	    //printf("saurabh\n");
	    double **temp2 = matrixMul(temp1, Inv2, x2-xmid+1);
	    Inv3 = negMatrix(temp2, x2-xmid+1, x2-xmid+1);

	    for(int i=0; i < xmid-x1; i++){
	        for(int j=0; j < ymid-y1; j++){
		        Inv[i][j] = Inv1[i][j];
		    }
	    }
	    for(int i=xmid-x1; i <= x2-x1; i++){
	        for(int j=ymid-y1; j<= y2-y1; j++){
		        Inv[i][j] = Inv2[i-xmid+x1][j-ymid+y1];
		    }
		}
		for(int i=0; i < xmid-x1; i++){
	        for(int j=ymid-y1; j<= y2-y1; j++){
		        Inv[i][j] = Inv3[i][j-ymid+y1];
		    }
		}
		for(int i=xmid-x1; i <= x2-x1; i++){
	        for(int j=ymid-y1; j<= y2-y1; j++){
		        Inv[i][j] = 0;
		    }
		}
		//printf("x1=%d y1=%d x2=%d y2=%d\n", x1, y1, x2, y2);
		//printMatrix(Inv, x2-x1+1, y2-y1+1);

		return(Inv);
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
    double** Inv = compute_inverse(Mat, 0, 0, n-1, n-1);
    //printf("compute_inverse time = %10.4e\n", omp_get_wtime() - start);
    //printMatrix(Mat, n, n, fp);
    printMatrix(Inv, n, n, fp);

    printf("\ncomputation norm = %lf\n", norm(Mat, Inv, n));
    fclose(fp);
    return(0);
}


