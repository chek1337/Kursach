#include "GRID.h"

int main()
{
	FILE* inFE, * inXY, * inZ, * inFirstBC, * inSecondBC, * inThirdBC;
	fopen_s(&inFE, "FE2.txt", "r");
	fopen_s(&inXY, "XY2.txt", "r");
	fopen_s(&inZ, "Z2.txt", "r");
	fopen_s(&inFirstBC, "inFirstBC2.txt", "r");
	fopen_s(&inSecondBC, "inSecondBC2.txt", "r");
	fopen_s(&inThirdBC, "inThirdBC2.txt", "r");
	GridAndSLAE grid;
	grid.InputFromFile(inFE, inXY, inZ, inFirstBC, inSecondBC, inThirdBC);

	grid.CalculateA_b();
	grid.SecondBoundaryConditions();
	grid.ThirdBoundaryConditions();
	grid.FirstBoundaryConditions();
	grid.OutputDense();

	grid.MSGForNonSymMatrixWithLuSqP();
	grid.OutputDense();
}