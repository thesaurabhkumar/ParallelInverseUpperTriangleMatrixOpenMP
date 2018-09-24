#include <iostream>
#include <stdlib.h>
#include <math.h>

void printMatrix(float** Mat, int x, int y)
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

float** allocMatrix(int x, int y){
    float **matrix = NULL;
    matrix = (float**)malloc(x * sizeof(float*));
    for(int i=0; i<x; i++){
        matrix[i] = (float*)malloc(y * sizeof(float));
    }
    return(matrix);
}

float** matrixMul(float** Mat1, float** Mat2, int n)
{
    float** res = allocMatrix(n, n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            //res[i][j] = 0;
            for(int k=0; k<n; k++){
                //printf("res=%f  %f  %f\n", res[i][j], Mat1[i][k], Mat2[k][j]);
                res[i][j] += Mat1[i][k]*Mat2[k][j];
                //printf("i=%d j=%d k=%d res=%f\n", i, j, k, res[i][j]);
            }
        }
    }
    //printMatrix(Mat1, n, n);
    //printMatrix(Mat2, n, n);
    //printMatrix(res, n, n);
    return(res);
}

float** seqInv(float **Mat, int x1, int y1, int x2, int y2){
    float **Inv = allocMatrix(x2 - x1 + 1, y2 - y1 + 1);
    float det = Mat[x1][y1]*Mat[x2][y2] - Mat[x1][y2]*Mat[x2][y1];
    Inv[0][0] = Mat[x2][y2]/det;
    Inv[0][1] = -Mat[x1][y2]/det;
    Inv[1][0] = -Mat[x2][y1]/det;
    Inv[1][1] = Mat[x1][y1]/det;
    //printf("x1=%d y1=%d x2=%d y2=%d\n", x1, y1, x2, y2);
	//printMatrix(Inv, x2-x1+1, y2-y1+1);
    return(Inv);
}

float** sliceMatrix(float** Mat, int x1, int y1, int x2, int y2)
{
    float** slice = allocMatrix(x2-x1+1, y2-y1+1);
    for(int i=0; i<x2-x1+1; i++){
        for(int j=0; j<y2-y1+1; j++){
            slice[i][j] = Mat[x1+i][y1+j];
        }
    }
    //printf("SLICE x1=%d y1=%d x2=%d y2=%d\n", x1, y1, x2, y2);
	//printMatrix(slice, x2-x1+1, y2-y1+1);
    return(slice);
}

float** negMatrix(float** Mat, int x, int y){
    float** res = allocMatrix(x+1, y+1);
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            res[i][j] = -Mat[i][j];
        }
    }
    return(res);
}

float** compute_inverse(float** Mat, int x1, int y1, int x2, int y2)
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
        float **Inv1 = NULL;
        float **Inv2 = NULL;
        float **Inv = allocMatrix(x2 - x1 + 1, y2 - y1 + 1);
        float **Inv3 = NULL;
        int xmid = round((x2 + x1 + 1)/2);
        int ymid = round((y2 + y1 + 1)/2);

	    Inv1 = compute_inverse(Mat, x1, y1, xmid-1, ymid-1);
	    Inv2 = compute_inverse(Mat, xmid, ymid, x2, y2);
        float** slice = sliceMatrix(Mat, x1, ymid, xmid-1, y2);
        //printf("RECUR x1=%d y1=%d xmid=%d ymid=%d x2=%d y2=%d\n", x1, y1, xmid, ymid, x2, y2);
        //printMatrix(Inv1, xmid-x1, ymid-y1);
        //printMatrix(Inv2, x2-xmid+1, y2-ymid+1);

        //printf("saurabh\n");
	    float **temp1 = matrixMul(Inv1, slice, x2-xmid+1);
	    //printf("saurabh\n");
	    float **temp2 = matrixMul(temp1, Inv2, x2-xmid+1);
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
		//printf("x1=%d y1=%d x2=%d y2=%d\n", x1, y1, x2, y2);
		//printMatrix(Inv, x2-x1+1, y2-y1+1);

		return(Inv);
	}	    

}

float** upperTrian(float** Mat, int n){
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

void initMatrix(float** Mat, int x, int y)
{
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            Mat[i][j] = rand()%20 + 1;
        }
    }
}

int main(int argc, char** argv){
    int n = 0;
    if(argc != 2){
        printf("Incorrect number of arguments: %d. Please Check!\n", argc);
        exit(-1);
    }
    n = atoi(argv[1]);
    float** Mat = allocMatrix(n, n);
    initMatrix(Mat, n, n);
    Mat = upperTrian(Mat,n);
    printMatrix(Mat, n, n);

    float** Inv = compute_inverse(Mat, 0, 0, n-1, n-1);
    printMatrix(Inv, n, n);


    return(0);
}


