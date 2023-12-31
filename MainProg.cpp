//#include "GRID.h"
#include "TESTS\TEST_aprox_x3_y3_z3\TEST_aprox_x3_y3_z3.h"
#include "string"

string path = "TESTS\\TEST_aprox_x3_y3_z3\\h2_XYZ\\";
string FE = path + "FE.txt";
string XY = path + "XY.txt";
string Z = path + "Z.txt";
string FBC = path + "FirstBC.txt";
string SBC = path + "SecondBC.txt";
string TBC = path + "ThirdBC.txt";

int main()
{
	FILE* inFE, * inXY, * inZ, * inFirstBC, * inSecondBC, * inThirdBC, *outQ;
	fopen_s(&inFE, FE.c_str(), "r");
	fopen_s(&inXY, XY.c_str(), "r");
	fopen_s(&inZ, Z.c_str(), "r");
	fopen_s(&inFirstBC, FBC.c_str(), "r");
	fopen_s(&inSecondBC, SBC.c_str(), "r");
	fopen_s(&inThirdBC, TBC.c_str(), "r");
	fopen_s(&outQ, "outQ.txt", "w");
	GridAndSLAE grid;
	grid.InputFromFile(inFE, inXY, inZ, inFirstBC, inSecondBC, inThirdBC);

	grid.CalculateA_b();
	grid.SecondBoundaryConditions();
	grid.ThirdBoundaryConditions();
	grid.FirstBoundaryConditions();
	//grid.OutputDense();

	grid.MSGForNonSymMatrixWithLuSqP();
	//grid.OutputDense();
	grid.OutputSolutionQ();
	//grid.OutputSolutionQinFile(outQ);

}