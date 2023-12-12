#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

#define REALOUTD "%.7f\t"
#define REALOUTDb "%.7f\n"
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
};

struct ThirdBoundary
{
    int node1, node2, node3, node4; // x1y1, x2y1, x1y2, x2y2
    int num_ubeta;
    double beta;
};

class GridAndSLAE
{
public:
    // NumberOfNodes = NoN
    int NoN; // Общее кол-во узлов
    int NoN_xy; // Общее кол-во узлов в плоскости XY
    int NoN_z; // Общее кол-во узлов на оси Z
    int NoN_fe; // Общее кол-во конечных элементов на оси Z
    int Nof_FisrtBC;
    int Nof_SecondBC;
    int Nof_ThirdBC;

    vector<FiniteElement> fe; // Массив содержащий информацию об о конечном элементе
    vector<Point> xy; // Массив с координатами сетки в Oxy
    vector<double> z; // Массив с координатами сетки вдоль оси Z
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
    
    void MSGForSymMatrixWithLuSqP();

protected:
    int maxiter = 10000;
    vector<int> ia;
    vector<int> ja;
    vector<double> al;
    vector<double> di;
    vector<double> x;
    vector<double> b;
    vector<double> r;
    vector<double> z_msg;
    vector<double> tmp;
    vector<double> x0;
    vector<double> alLU;
    vector<double> diLU;
    vector<vector<int>> iaja;
    double eps = 1e-13;
    int maxIter = 10000;


    void GeneratePortrait();
    void VectorCopy(vector<double>& from, vector<double>& to);
    void CalculateRelativeDiscrepancy(vector<double>& vectorMult, vector<double>& vectorOut);
    void MatrixVectorMultiplication(vector<double>& vectorMult, vector<double>& vectorOut);
    void SolveForwardLU(vector<double>& lowerTringMat, vector<double>& diag, vector<double>& rightVector, vector<double>& vectorX);
    void SolveBackwardLU(vector<double>& upperTringMat, vector<double>& diag, vector<double>& rightVector, vector<double>& vectorX);
    void CalculateLUsq();
    double CalculateRelativeDiscrepancyWithR(double norm);
    double VectorNorm(vector<double>& vector);
    double VectorScalarProduction(vector<double>& vector1, vector<double>& vector2);

    double Gauss3_Gxy(int i, int j, double b1, double b2, double b3, double b4, double b5, double b6, double a0, double a1, double a2);
    double Phi(double e, double n, int i, int j, double b1, double b2, double b3, double b4, double b5, double b6, double a0, double a1, double a2);
    double TETA(int number, double x, double y, double z);
    double UBETA(int number, double x, double y, double z);
    void AllocateMemory();
};

