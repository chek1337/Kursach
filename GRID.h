#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

#define REALOUTD "%.4lf\t"
#define REALOUTDb "%.4lf\n"
using namespace std;

struct Point {
    double x;
    double y;
};

struct FiniteElement
{
    int node1, node2, node3, node4; // x1y1, x2y1, x1y2, x2y2
    int bottom, top; //z1, z2
    int region;
};

struct FirstBoundary
{
    int node;
    double ug;
};

struct SecondBoundary
{
    int node1, node2, node3, node4; // x1y1, x2y1, x1y2, x2y2
    int num_teta;
    int side; // -1 - ���������, 1 - ������� �����
};

struct ThirdBoundary
{
    int node1, node2, node3, node4; // x1y1, x2y1, x1y2, x2y2
    int num_ubeta;
    double beta;
    int side; // 0 - ���������, 1 - ������� �����
};

class GridAndSLAE
{
public:
    // NumberOfNodes = NoN
    int NoN; // ����� ���-�� �����
    int NoN_xy; // ����� ���-�� ����� � ��������� XY
    int NoN_z; // ����� ���-�� ����� �� ��� Z
    int NoN_fe; // ����� ���-�� �������� ��������� �� ��� Z
    int Nof_FisrtBC;
    int Nof_SecondBC;
    int Nof_ThirdBC;

    vector<FiniteElement> fe; // ������ ���������� ���������� �� � �������� ��������
    vector<Point> xy; // ������ � ������������ ����� � Oxy
    vector<double> z; // ������ � ������������ ����� ����� ��� Z
    vector<FirstBoundary> FirstBC;
    vector<SecondBoundary> SecondBC;
    vector<ThirdBoundary> ThirdBC;

    void InputFromFile(FILE* inFE, FILE* inXY, FILE* inZ, FILE* inFirstBC, FILE* inSecondBC, FILE* inThirdBC);
    void CalculateA_b();
    void FirstBoundaryConditions();
    void SecondBoundaryConditions();
    void ThirdBoundaryConditions();
    void OutputDense();
    void OutputLUDense();
    //void SolveSLAE();
    void MatrixVectorMultiplication(vector<double>& vectorMult, vector<double>& vectorOut);
    void MSGForNonSymMatrixWithLuSqP();
    void OutputSolutionQ();
    vector<double> x;
    vector<double> tmp;
protected:
    int maxiter = 10000;
    vector<int> ia;
    vector<int> ja;
    vector<double> al;
    vector<double> au;
    vector<double> di;
   
    vector<double> b;
    vector<double> r;
    vector<double> z_msg;
   
    vector<double> x0;
    vector<double> alLU;
    vector<double> auLU;
    vector<double> diLU;
    vector<vector<int>> iaja;
    double eps = 1e-20;
    int maxIter = 10000;


    void GeneratePortrait();
    void VectorCopy(vector<double>& from, vector<double>& to);
    void CalculateRelativeDiscrepancy(vector<double>& vectorMult, vector<double>& vectorOut);
    
    void SolveForwardLU(vector<double>& lowerTringMat, vector<double>& diag, vector<double>& rightVector, vector<double>& vectorX);
    void SolveBackwardLU(vector<double>& upperTringMat, vector<double>& diag, vector<double>& rightVector, vector<double>& vectorX);
    void CalculateLUsq();
    double CalculateRelativeDiscrepancyWithR(double norm);
    double VectorNorm(vector<double>& vector);
    double VectorScalarProduction(vector<double>& vector1, vector<double>& vector2);

   
    void MatrixUVectorMultiplicationLU(vector<double>& upperTringMat, vector<double>& diag, vector<double>& vectorMult, vector<double>& vectorOut);
    void TransposedMatrixVectorMultiplication(vector<double>& vectorMult, vector<double>& vectorOut);
    void CalculateZ_LUsq(vector<double>& vectorOut);
    double CalculateRelativeDiscrepancy(double norm);
    void VectorSubtract(vector<double>& first, vector<double>& second, vector<double>& result);


    double Gauss3_Gxy(int i, int j, double b1, double b2, double b3, double b4, double b5, double b6, double a0, double a1, double a2);
    double Phi(double e, double n, int i, int j, double b1, double b2, double b3, double b4, double b5, double b6, double a0, double a1, double a2);
    double TETA(int number, double x, double y, double z);
    double UBETA(int number, double x, double y, double z);
    double GAMMA(int number);
    double LAMBDA(int number);
    double FUN(int number, double x, double y, double z);
    void AllocateMemory();
};

