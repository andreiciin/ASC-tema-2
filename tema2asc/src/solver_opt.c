/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	
	// C = B x A x At + Bt x B -- operatia
	// aloc memorie pt matricile in care salvez etapele operatiei
	register double *BxA = calloc(N * N, sizeof(double));
	register double *BtxB = calloc(N * N, sizeof(double));
	register double *BxAxAt = calloc(N * N, sizeof(double));
	// matricea finala pe care o returnez
	register double *C = calloc(N * N, sizeof(double));

	// Bt x B
        for (register int i = 0; i != N; i++) {
		// asociez pt fiecare matrice din adunarea:
		// BtxB[i * N + j] += B[k * N + j] * B[k * N + i] un pointer
                register double *op2_p = BtxB;
		// localizez linia i in matrice
                op2_p += i * N;
                register double *og_pB = B;

		// coloana din BtxB e incrementata in forul urmator
                for (register int j = 0; j != N; j++, op2_p += 1) {
                        register double mu = 0;

                        register double *pBt = og_pB;
			// localizez coloana j din B transpus
			pBt += j;
                        register double *pB = B;
			// localizez coloana i din B
			pB += i;
			// liniile lui B si Bt sunt incrementate in for
                        for (register int k = 0; k != N; k++, pBt += N, pB += N) {
                                mu += *pBt * *pB;
                        }
                        *op2_p += mu;
                }
        }
		
	// B x A
	for (register int i = 0; i != N; i++) {
		// cream pointeri pt urmatoarele matrice:
		// BxA[i * N + j] += B[i * N + k] * A[k * N + j]
		// localizez liniile prin i
		register double *og_p = BxA;
		og_p += i * N;
		register double *og_pB = B;
		og_pB += i * N;

		// coloana din BxA e incrementata in acest for
		for (register int j = 0; j != N; j++, og_p += 1) {
			register double mu = 0;
			
			register double *pB = og_pB;
			register double *pA = A;
			// localizez coloana j din A
			pA += j;
			
			// coloana din B si linia din A sunt actualizate in for
			for (register int k = 0; k != j + 1; k++, pB += 1, pA += N) {
				mu += *pB * *pA;
			}
			*og_p += mu;
		}
	}

	// B x A x At
	for (register int i = 0; i != N; i++) {
		// cream pointeri pt urmatoarele matrice:
		// BxAxAt[i * N + j] += BxA[i * N + k] * A[j * N + k];
		// localizez liniile prin i
		register double *op1_p = BxAxAt;
		op1_p += i * N;
		register double *og_pLeft = BxA;
		og_pLeft += i * N;

		// coloana lui BxAxAt e incrementata in for
		for (register int j = 0; j != N; j++, op1_p += 1) {
			register double mu = 0;
			
			register double *pLeft = og_pLeft;
			// incrementez coloana cu j deoarece k incepe de la j
			pLeft += j;
			register double *pAt = A;
			// linia j + coloana de la j
			pAt += N * j + j;
			
			for (register int k = j; k != N; k++, pLeft += 1, pAt += 1) {
				mu += *pLeft * *pAt;
			}
			*op1_p += mu;
		}
	}			

	// C = B x A x At + Bt x B
	for (register int i = 0; i < N; i++)
		for (register int j = 0; j < N; j++)
			*(C + i * N + j) += *(BxAxAt + i * N + j) + *(BtxB + i * N + j);

	free(BxA);
	free(BtxB);
	free(BxAxAt);
	return C;
}
