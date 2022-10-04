/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {

	// C = B x A x At + Bt x B  -- operatia
	// aloc memorie pt matricile in care salvez etapele operatiei
	double *BxA = calloc(N * N, sizeof(double));
	double *BtxB = calloc(N * N, sizeof(double));
	double *BxAxAt = calloc(N * N, sizeof(double));
	// matricea finala pe care o returnez
	double *C = calloc(N * N, sizeof(double));
	
	// adunarea e comutativa
	// Bt x B
        for (int i = 0; i != N; ++i)
                for (int j = 0; j != N; ++j)
                        for (int k = 0; k != N; ++k)
                                *(BtxB + i * N + j) += *(B + k * N + j) * *(B + k * N + i);
	// B x A
	for (int i = 0; i != N; i++)
		for (int j = 0; j != N; j++)
			for (int k = 0; k != j + 1; k++)
				*(BxA + i * N + j) += *(B + i * N + k) * *(A + k * N + j);

	// B x A x At
	for (int i = 0; i != N; i++)
		for (int j = 0; j != N; j++)
			for (int k = j; k != N; k++)
				*(BxAxAt + i * N + j) += *(BxA + i * N + k) * *(A + j * N + k);

	// C = B x A x At + Bt x B
	for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
				*(C + i * N + j) += *(BxAxAt + i * N + j) + *(BtxB + i * N + j); 
	
	free(BxA);
	free(BtxB);
	free(BxAxAt);
	return C;
}
