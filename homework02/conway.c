#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int par_num;
int iteration;
int xNum;
int yNum;
int dotsNum;

void initGame(matrix,x,y)
	int **matrix;
	int *x; 
	int *y;
{
	printf("initializing game...\n");
	int i , j;
	int tmpx, tmpy;
	for (i=0; i < yNum; i++){
		 for (j=0; j < xNum; j++){
			*((int*)matrix+i*xNum+j) = 0;
		}
	}
	printf("map initialization done!\n");
	for( i =0; i< dotsNum;i++){
		tmpy = *y;
		tmpx = *x;
		*((int*)matrix+ tmpy * xNum + tmpx) = 1;
		x++;
		y++;
	}
}

int countNeighbors(submatrix, xNum, yNum, x, y )
	int *submatrix;
	int xNum;
	int yNum;
	int x;
	int y;
{
	int count =0 ;

		if(   x>0 && y>0 && *(submatrix + (y-1)*xNum + (x-1)) == 1 ){
			count++;
		}
		if(   y>0 && *(submatrix + (y-1)*xNum + (x)) == 1 ){
			count++;
		}
		if(  y>0 && x<(xNum-1 ) && *(submatrix + (y-1)*xNum + (x+1)) == 1  ){
			count++;
		}
		if( x>0 && *(submatrix + (y)*xNum + (x-1)) == 1  ){
			count++;
		}
		if( x<(xNum-1) && *(submatrix + (y)*xNum + (x+1)) == 1 ){
			count++;
		}
		if( y<(yNum-1) && x>0 &&  *(submatrix + (y+1)*xNum + (x-1)) == 1 ){
			count++;
		}
		if(  y<(yNum-1) && *(submatrix + (y+1)*xNum + (x)) == 1 ){
			count++;
		}
		if( y<(yNum-1) && x<(xNum-1) && *(submatrix + (y+1)*xNum + (x+1)) == 1){
			count++;
		}
	
	return count;
			 

}

void printMatrix(matrix, xNum, yNum)
	int **matrix;
	int xNum;
	int yNum;
{
	int i,j;
	for(i=0; i<yNum; i++){		
		for(j=0; j<xNum; j++){
			printf(" %d ", *((int *)matrix+ i* xNum + j));
		}
		printf("\n");
	}
}

int main (argc, argv)
	int argc;
	char *argv[];  //usage: <par-NUM> <xNUM> <yNUM> <dotsNUM> <x1> <y1>....<xn> <yn>
{
	if(argc < 4){
		printf("Error Parameters\n");
		printf(" //usage: ./conway <par-NUM> <Iteration> <xNUM> <yNUM> <dotsNUM> <x1> <y1>....<xn> <yn>\n");
		return -1;
	}

	int rank, numtasks, rc;

	rc = MPI_Init(&argc, &argv);
	if(rc != MPI_SUCCESS){
		printf("Error: cannot start mpi program!\n");
		MPI_Abort(MPI_COMM_WORLD,rc);
	} 
	MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
		
	
	par_num = atoi(argv[1]);
	iteration = atoi(argv[2]);
	xNum = atoi(argv[3]);
	yNum = atoi(argv[4]);
	int matrix[yNum][xNum];
	dotsNum = atoi(argv[5]);
	int i,j ;
	int x[dotsNum], y[dotsNum];
	for ( i =0 ; i<dotsNum; i++){
		x[i] = atoi(argv[6+2*i]);
		y[i] = atoi(argv[6+2*i+1]);
	}
	
	initGame(matrix, x, y);
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0 )
	{
		printMatrix(matrix, xNum, yNum);
	}

		
	if(numtasks != par_num){
		printf(" Error:  Parallel tasks !\n");
		printf(" Num of available tasks: %d \n", numtasks);
		return -1;
	}
	if(par_num < 2){
		printf(" Error: too less parallel tasks !\n");
		printf(" Please set Par-Num > 1\n");
		return -1;
	}
	
	
	int time = 1;
	int *subMatrix = (int *)matrix;
	int subYNum = yNum / par_num;
	int subYStart = subYNum * rank;
	int subYEnd = subYNum * (rank+1) - 1;
	if(rank == par_num -1)
	{
		subYEnd = yNum-1;
		subYNum = yNum - subYStart ;
	}
	int tmpMatrix[subYNum][xNum];	
	int upBoundSendBuf[xNum];
	int downBoundSendBuf[xNum];
	int upBoundRecvBuf[xNum];
	int downBoundRecvBuf[xNum]; 
	int tagMsgUpToDown =1;
	int tagMsgDownToUp =2;
	int count = 0;
	for(time = 1; time <= iteration; time++){
		MPI_Request reqs01[2];
		MPI_Request reqs02[2];
		MPI_Status stats01[2];
		MPI_Status stats02[2];
		if(rank == 0){
			printf("iteration = %d ..\n", time);
		}
		if(rank != (par_num-1)){
			MPI_Irecv(&downBoundRecvBuf, xNum, MPI_INT, rank+1, tagMsgDownToUp, MPI_COMM_WORLD, &reqs01[0]);
			for(i = 0; i < xNum; i++){
				downBoundSendBuf[i] = *( subMatrix + subYEnd * xNum + i);
			}
			MPI_Isend(&downBoundSendBuf, xNum, MPI_INT, rank+1, tagMsgUpToDown, MPI_COMM_WORLD, &reqs01[1]);	
		}
		if(rank != 0 ){
			MPI_Irecv(&upBoundRecvBuf, xNum, MPI_INT, rank-1, tagMsgUpToDown, MPI_COMM_WORLD, &reqs02[0]);
			for(i = 0; i< xNum; i++){
				upBoundSendBuf[i] = *( subMatrix + subYStart * xNum + i);
			}
			MPI_Isend(&upBoundSendBuf, xNum, MPI_INT, rank-1, tagMsgDownToUp, MPI_COMM_WORLD, &reqs02[1]);
		}
		
		//start count the neighbors for the inner dots
		int tmpBegin, tmpEnd;
		if(rank == 0) 
		{
			tmpBegin = subYStart;
			tmpEnd = subYEnd - 1;
		}
		if(rank == par_num -1 )
		{
			tmpBegin = subYStart + 1;
			tmpEnd = subYEnd;
		}
		if( rank > 0 && rank < par_num-1 )
		{
			tmpBegin = subYStart + 1;
			tmpEnd = subYEnd - 1;
		}
		for(i = tmpBegin ; i <= tmpEnd ; i++){
			for( j = 0 ; j < xNum ; j++ ){
				count =	countNeighbors(subMatrix, xNum, yNum, j ,i );				
				if(count < 2){
					tmpMatrix[i-subYStart][j] = 0;
				}
				if(count== 2 ) {
					if( *(subMatrix+i*xNum+j)==1){
						tmpMatrix[i-subYStart][j] = 1;
					}
					else {
						tmpMatrix[i-subYStart][j] = 0;
					}
				}
				if(count == 3 ){
					tmpMatrix[i-subYStart][j] = 1;
				}
				if(count>3) {
					tmpMatrix[i-subYStart][j] = 0;
				}
			}
		}

		
		if(rank != par_num - 1){
			MPI_Waitall(2, reqs01, stats01 );
		}
		if(rank != 0){
			MPI_Waitall(2, reqs02, stats02 );
		}

	
		//after communication done. calculate the num of live neighbors for the bound
		if( rank != (par_num-1) ){
			for ( j = 0; j < xNum; j++){
				*( subMatrix + (subYEnd+1)*xNum + j) = downBoundRecvBuf[j];	
			}
			i = subYEnd;
			for ( j = 0; j < xNum; j++){
				count = countNeighbors(subMatrix, xNum, yNum, j, i);
				if(count < 2){
					tmpMatrix[i-subYStart][j] = 0;
				}
				if(count== 2 ) {
					if( *(subMatrix+i*xNum+j)==1){
						tmpMatrix[i-subYStart][j] = 1;
					}
					else {
						tmpMatrix[i-subYStart][j] = 0;
					}
				}
				if(count == 3 ){
					tmpMatrix[i-subYStart][j] =1;
				}
				if(count>3) {
					tmpMatrix[i-subYStart][j] = 0;
				}
			}
		}

		if( rank != 0 ){
			for( j = 0; j < xNum; j++){
				*( subMatrix + (subYStart-1)*xNum +j) = upBoundRecvBuf[j];
			}
			i = subYStart;
			for( j = 0; j < xNum; j++){
				count = countNeighbors(subMatrix, xNum, yNum, j, i);
				if(count < 2){
					tmpMatrix[i-subYStart][j] = 0;
				}
				if(count== 2 ) {
					if( *(subMatrix+i*xNum+j)==1){
						tmpMatrix[i-subYStart][j] = 1;
					}
					else {
						tmpMatrix[i-subYStart][j] = 0;
					}
				}
				if(count == 3 ){
					tmpMatrix[i-subYStart][j] =1;
				}
				if(count>3) {
					tmpMatrix[i-subYStart][j] = 0;
				}
			}
		}

		for(i = subYStart; i <= subYEnd; i++){
			for(j = 0; j < xNum; j++){
				*(subMatrix + i*xNum + j) = tmpMatrix[i-subYStart][j];
			}
		}

		//sequencial print for each iteration
		int k;
		char tmpC;
		for( k=0; k< par_num;k++){
			if(rank == k){		
				for(i = subYStart; i <= subYEnd; i++){
					for(j = 0; j < xNum; j++){
						if( *(subMatrix + i*xNum +j) == 1 ){
							tmpC = 'X';
						}
						else tmpC = '-';
						printf(" %c ", tmpC);
					}
					printf("	line: %d\n",i);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}
	
	MPI_Finalize();	 
	return 0;
}
