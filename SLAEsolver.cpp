#include "GRID.h"

void GridAndSLAE::OutputDense()
{
	int flagfound = 0;
	for (int i = 0; i < NoN; i++)
	{
		int k = ia[i + 1] - ia[i];
		if (k == 0)
		{
			for (int j = 0; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}
		else
		{
			int lastj = 0;
			for (int j = ia[i]; j < ia[i + 1]; j++)
			{
				for (int p = lastj; p < ja[j]; p++)
				{
					printf(REALOUTD, 0.0);
				}
				printf(REALOUTD, al[j]);
				lastj = ja[j] + 1;
			}
			for (int j = lastj; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}

		printf(REALOUTD, di[i]);

		for (int j = i + 1; j < NoN; j++)
		{
			k = ia[j + 1] - ia[j];
			if (k == 0) {
				printf(REALOUTD, 0.0);
			}
			else
			{
				flagfound = 0;
				for (k = ia[j]; k < ia[j + 1]; k++)
				{

					if (ja[k] == i)
					{
						printf(REALOUTD, al[k]);
						flagfound = 1;
						break;
					}
				}
				if (flagfound == 0)
					printf(REALOUTD, 0.0);
			}
		}
		printf("\n");
	}
	printf("\n");

	for (int i = 0; i < NoN; i++)
	{
		printf(REALOUTDb, b[i]);
	}
}

void GridAndSLAE::OutputLUDense()
{
	int flagfound = 0;
	for (int i = 0; i < NoN; i++)
	{
		int k = ia[i + 1] - ia[i];
		if (k == 0)
		{
			for (int j = 0; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}
		else
		{
			int lastj = 0;
			for (int j = ia[i]; j < ia[i + 1]; j++)
			{
				for (int p = lastj; p < ja[j]; p++)
				{
					printf(REALOUTD, 0.0);
				}
				printf(REALOUTD, alLU[j]);
				lastj = ja[j] + 1;
			}
			for (int j = lastj; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}

		printf(REALOUTD, diLU[i]);

		for (int j = i + 1; j < NoN; j++)
		{
			k = ia[j + 1] - ia[j];
			if (k == 0) {
				printf(REALOUTD, 0.0);
			}
			else
			{
				flagfound = 0;
				for (k = ia[j]; k < ia[j + 1]; k++)
				{

					if (ja[k] == i)
					{
						printf(REALOUTD, alLU[k]);
						flagfound = 1;
						break;
					}
				}
				if (flagfound == 0)
					printf(REALOUTD, 0.0);
			}
		}
		printf("\n");
	}
	printf("\n");
}


void GridAndSLAE::VectorCopy(vector<double>& from, vector<double>& to)
{
	for (int i = 0; i < NoN; i++)
		to[i] = from[i];
}

void GridAndSLAE::CalculateRelativeDiscrepancy(vector<double>& vectorMult, vector<double>& vectorOut)
{
	for (int i = 0; i < NoN; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		vectorOut[i] = b[i] - vectorOut[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			vectorOut[i] += al[k] * vectorMult[j];
			vectorOut[j] += al[k] * vectorMult[i];
		}
	}
}


void GridAndSLAE::MatrixVectorMultiplication(vector<double>& vectorMult, vector<double>& vectorOut)
{
	for (int i = 0; i < NoN; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			vectorOut[i] += al[k] * vectorMult[j];
			vectorOut[j] += al[k] * vectorMult[i];
		}
	}
}

void GridAndSLAE::SolveForwardLU(vector<double> &lowerTringMat, vector<double> &diag, vector<double> &rightVector, vector<double> &vectorX) {
	for (int i0, i1, i = 0; i < NoN; i++)
	{
		double sum = 0;
		i0 = ia[i];
		i1 = ia[i + 1];
		//int j = i - (i1 - i0);
		for (int k = i0; k < i1; k++)
		{
			int j = ja[k];
			sum += lowerTringMat[k] * vectorX[j];
		}
		vectorX[i] = (rightVector[i] - sum) / diag[i];
	}
}

void GridAndSLAE::SolveBackwardLU(vector<double>& upperTringMat, vector<double>& diag, vector<double>& rightVector, vector<double>& vectorX) {
	for (int i0, i1, i = NoN - 1; i >= 0; i--)
	{
		i0 = ia[i];
		i1 = ia[i + 1];
		vectorX[i] = rightVector[i] / diag[i];
		for (int j, k = i0; k < i1; k++)
		{
			j = ja[k];
			rightVector[j] -= upperTringMat[k] * vectorX[i];
		}
	}
}

double GridAndSLAE::VectorNorm(vector<double>& vector)
{
	double norm = 0;
	for (int i = 0; i < NoN; i++)
		norm += vector[i] * vector[i];
	return sqrt(norm);
}

double GridAndSLAE::VectorScalarProduction(vector<double>& vector1, vector<double>& vector2)
{
	double prod = 0;
	for (int i = 0; i < NoN; i++)
	{
		prod += vector1[i] * vector2[i];
	}
	return prod;
}

double GridAndSLAE::CalculateRelativeDiscrepancyWithR(double norm)
{
	return VectorNorm(r) / norm;
}

void GridAndSLAE::CalculateLUsq()
{
	for (int i = 0; i < NoN; i++) //i = 7
	{
		int i0 = ia[i]; //15
		int i1 = ia[i + 1]; //19
		double sumD = 0;

		for (int k = i0; k < i1; k++)
		{
			double sumL = 0, sumU = 0;
			int j = ja[k];

			int j0 = ia[j];
			int j1 = ia[j + 1];

			int ik = i0;
			int kj = j0;

			for (; ik < k and kj < j1; )
			{
				if (ja[ik] < ja[kj])
					ik++;
				if (ja[ik] > ja[kj])
					kj++;
				if (ja[ik] == ja[kj]) {
					sumL += alLU[ik] * alLU[kj];
					ik++; kj++;
				}

			}
			alLU[k] = (al[k] - sumL) / diLU[j];
			sumD += alLU[k] * alLU[k];
		}
		diLU[i] = sqrt(di[i] - sumD);
	}
}


void GridAndSLAE::MSGForSymMatrixWithLuSqP()
{
	CalculateLUsq();
	//OutputLUDense();
	VectorCopy(x0, x);
	CalculateRelativeDiscrepancy(x0, r);
	// M = LU
	// M^(-1) = U^(-1)*L^(-1)
	SolveForwardLU(alLU, diLU, r, tmp);
	SolveBackwardLU(alLU, diLU, tmp, z);


	double normB = VectorNorm(b);
	double ak = 0, bk = 0;

	//ak
	SolveForwardLU(alLU, diLU, r, tmp);
	SolveBackwardLU(alLU, diLU, tmp, tmp);
	double Mr_rPrev = VectorScalarProduction(tmp, r); // ( M^(-1)*r(k-1) , r(k-1) )
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		MatrixVectorMultiplication(z, tmp); // A*z(k-1)
		double Az_zPrev = VectorScalarProduction(tmp, z); //( A*z(k-1) , z(k-1) )
		ak = Mr_rPrev / Az_zPrev;

		for (int i = 0; i < NoN; i++)
		{
			x[i] = x[i] + ak * z[i]; // xk = x(k-1) + ak*z(k-1)
			r[i] = r[i] - ak * tmp[i]; // rk = r(k-1) - ak*U^(-T)*A^T*L^(-T)*L^(-1)*A*U^(-1)*z(k-1)
		}

		SolveForwardLU(alLU, diLU, r, tmp);
		SolveBackwardLU(alLU, diLU, tmp, tmp);
		double Mr_rCur = VectorScalarProduction(tmp, r); // ( rk , rk )

		bk = Mr_rCur / Mr_rPrev;

		for (int i = 0; i < NoN; i++) // zk =  M^(-1)*rk + bk*z(k-1)
			z[i] = tmp[i] + bk * z[i];

		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		Mr_rPrev = Mr_rCur;
		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}
}

