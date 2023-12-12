#include "GRID.h"

int main()
{
	FILE* inFE, * inXY, * inZ, * inFirstBC, * inSecondBC, * inThirdBC;
	fopen_s(&inFE, "FE1.txt", "r");
	fopen_s(&inXY, "XY1.txt", "r");
	fopen_s(&inZ, "Z1.txt", "r");
	fopen_s(&inFirstBC, "inFirstBC1.txt", "r");
	fopen_s(&inSecondBC, "inSecondBC.txt", "r");
	fopen_s(&inThirdBC, "inThirdBC.txt", "r");
	GridAndSLAE grid;
	grid.InputFromFile(inFE, inXY, inZ, inFirstBC, inSecondBC, inThirdBC);

	grid.CalculateA_b();
	grid.SecondBoundaryConditions();
	grid.ThirdBoundaryConditions();
	grid.OutputDense();

	grid.FirstBoundaryConditions();
	grid.OutputDense();

	grid.MSGForSymMatrixWithLuSqP();
	grid.OutputDense();
}