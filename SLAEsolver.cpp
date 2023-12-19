#include "GRID.h"

//void GridAndSLAE::Input(FILE* paramf, FILE* iaf, FILE* jaf, FILE* alf, FILE* auf, FILE* dif, FILE* bf) {
//	fscanf_s(paramf, "%d", &NoN);
//	fscanf_s(paramf, "%d", &maxiter);
//	fscanf_s(paramf, "%lf", &eps);
//
//	ia.resize(NoN + 1);
//	for (int i = 0; i <= NoN; i++)
//		fscanf_s(iaf, "%d", &ia[i]);
//	int nProfile = ia[NoN] - ia[0];
//	ja.resize(nProfile);
//	for (int i = 0; i < nProfile; i++)
//		fscanf_s(jaf, "%d", &ja[i]);
//
//	if (ia[0]) {
//		for (int i = 0; i <= NoN; i++)
//			ia[i]--;
//		for (int i = 0; i < nProfile; i++)
//			ja[i]--;
//	}
//
//	al.resize(nProfile);
//	alLU.resize(nProfile);
//	for (int i = 0; i < nProfile; i++)
//		fscanf_s(alf, "%lf", &al[i]);
//
//	au.resize(nProfile);
//	auLU.resize(nProfile);
//	for (int i = 0; i < nProfile; i++)
//		fscanf_s(auf, "%lf", &au[i]);
//
//	di.resize(NoN);
//	diLU.resize(NoN);
//	for (int i = 0; i < NoN; i++)
//		fscanf_s(dif, "%lf", &di[i]);
//	b.resize(NoN);
//	for (int i = 0; i < NoN; i++)
//		fscanf_s(bf, "%lf", &b[i]);
//
//
//	x.resize(NoN);
//	r.resize(NoN);
//	z.resize(NoN);
//	tmp.resize(NoN);
//	x0.resize(NoN);
//	//for (int i = 0; i < n; i++)
//	//	x0[i] = 0;
//}
//


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
				printf(REALOUTD, al[j]); // ��� ������� 216
			    //printf(REALOUTD, 216 * al[j]); // ��� ������� 216
				lastj = ja[j] + 1;
			}
			for (int j = lastj; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}

		printf(REALOUTD, di[i]); // ��� ������� 216
		//printf(REALOUTD, 216 * di[i]); // ��� ������� 216

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
						printf(REALOUTD, au[k]); // ��� ������� 216
						//printf(REALOUTD, 216 * au[k]); // ��� ������� 216
						flagfound = 1;
						break;
					}
				}
				if (flagfound == 0)
					printf(REALOUTD, 0.0);
			}
		}
		printf("|   %.4lf\n", b[i]); // ��� ������� 216
		//printf("|   %.4lf\n", 216 * b[i]); // ��� ������� 216
		printf("\n");
	}
	printf("\n");
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
						printf(REALOUTD, auLU[k]);
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
			vectorOut[j] += au[k] * vectorMult[i];
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
			vectorOut[j] += au[k] * vectorMult[i];
		}
	}
}

void GridAndSLAE::SolveForwardLU(vector<double>& lowerTringMat, vector<double>& diag, vector<double>& rightVector, vector<double>& vectorX) {
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
					sumL += alLU[ik] * auLU[kj];
					sumU += alLU[kj] * auLU[ik];
					ik++; kj++;
				}

			}
			alLU[k] = (al[k] - sumL) / diLU[j];
			auLU[k] = (au[k] - sumU) / diLU[j];
			sumD += alLU[k] * auLU[k];
		}
		diLU[i] = sqrt(di[i] - sumD);
	}
}



void GridAndSLAE::MatrixUVectorMultiplicationLU(vector<double>& upperTringMat, vector<double>& diag, vector<double>& vectorMult, vector<double>& vectorOut)
{
	for (int i = 0; i < NoN; i++)
	{
		vectorOut[i] = vectorMult[i] * diag[i];
		for (int j, k = ia[i]; k < ia[i + 1]; k++)
		{
			j = ja[k];
			vectorOut[j] += upperTringMat[k] * vectorMult[i];
		}
	}
}

void GridAndSLAE::TransposedMatrixVectorMultiplication(vector<double>& vectorMult, vector<double>& vectorOut)
{
	for (int i = 0; i < NoN; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			vectorOut[i] += au[k] * vectorMult[j];
			vectorOut[j] += al[k] * vectorMult[i];
		}
	}
}

void GridAndSLAE::CalculateZ_LUsq(vector<double>& vectorOut)
{
	SolveBackwardLU(auLU, diLU, tmp, tmp);
	MatrixVectorMultiplication(tmp, x0);
	SolveForwardLU(alLU, diLU, x0, tmp);
	SolveBackwardLU(alLU, diLU, tmp, tmp);
	TransposedMatrixVectorMultiplication(tmp, x0);
	SolveForwardLU(auLU, diLU, x0, vectorOut);
}

void GridAndSLAE::VectorSubtract(vector<double>& first, vector<double>& second, vector<double>& result)
{
	for (int i = 0; i < NoN; i++)
	{
		result[i] = first[i] - second[i];
	}
}

double GridAndSLAE::CalculateRelativeDiscrepancy(double norm)
{
	MatrixVectorMultiplication(x, tmp);
	VectorSubtract(b, tmp, tmp);
	return VectorNorm(tmp) / norm;
}

void GridAndSLAE::OutputSolutionQ()
{
	printf("q = {\n");
	for (int i = 0; i < NoN - 1; i++)
	{
		printf("%.15lf\n", x[i]);
	}
	printf("%.15lf\n }\n", x[NoN - 1]);
}

void GridAndSLAE::MSGForNonSymMatrixWithLuSqP()
{
	CalculateLUsq();
	MatrixUVectorMultiplicationLU(auLU, diLU, x0, x);

	CalculateRelativeDiscrepancy(x0, tmp);
	SolveForwardLU(alLU, diLU, tmp, tmp);
	SolveBackwardLU(alLU, diLU, tmp, tmp);
	TransposedMatrixVectorMultiplication(tmp, x0);
	SolveForwardLU(auLU, diLU, x0, tmp);
	VectorCopy(tmp, r);
	VectorCopy(tmp, z);

	double r_rPrev, r_rCur, Newz_zPrev, ak, bk;
	r_rPrev = VectorScalarProduction(r, r);
	double normB = VectorNorm(b);
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		VectorCopy(z, tmp);
		CalculateZ_LUsq(tmp);
		Newz_zPrev = VectorScalarProduction(tmp, z);
		ak = r_rPrev / Newz_zPrev;

		for (int i = 0; i < NoN; i++)
		{
			x[i] = x[i] + ak * z[i];
			r[i] = r[i] - ak * tmp[i];
		}

		r_rCur = VectorScalarProduction(r, r);
		bk = r_rCur / r_rPrev;

		for (int i = 0; i < NoN; i++)
			z[i] = r[i] + bk * z[i];

		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		r_rPrev = r_rCur;

		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}
	SolveBackwardLU(auLU, diLU, x, x);
	printf("%.15lf\n", CalculateRelativeDiscrepancy(normB));

	MatrixVectorMultiplication(x, tmp);

}