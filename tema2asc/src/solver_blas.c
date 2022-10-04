/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <cblas.h>

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {

	// C = B x A x At + Bt x B -- operatia
	// aloc memorie pt matricile in care salvez etapele operatiei
	double *BtxB = calloc(N * N, sizeof(double));
	double *AxAt = calloc(N * N, sizeof(double));
	double *BxAxAt = calloc(N * N, sizeof(double));
	// matricea finala pe care o returnez
	double *C = calloc(N * N, sizeof(double));

	// Bt x B
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0,
                B, N, B, N, 1.0, BtxB, N);

	// tratez matricele triunghiulare
	// A x At
	for (int i = 0; i != N; i++)
		for (int j = 0; j != N; j++)
			*(AxAt + i * N + j) = *(A + i * N + j);	

	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit,
		N, N, 1.0, A, N, AxAt, N);	

	//  B x A x At
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0,
		B, N, AxAt, N, 0.0, BxAxAt, N);

	// B x A x At + Bt x B
	// am cautat si adunarea matricelor nu exista in rutinele blas
	for (int i = 0; i != N; i++)
		for (int j = 0; j != N; j++) 
			*(C + i * N + j) = *(BxAxAt + i * N + j) + *(BtxB + i * N + j);


	free(BtxB);
	free(AxAt);
	free(BxAxAt);

	return C;
	
}
