#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
using namespace std;

struct Point {
    double x;
    double y;
};

struct FiniteElement
{
    int node1, node2, node3, node4; // x1y1, x2y1, x1y2, x2y2
    int bottom, top; //z1, z2
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
    int num_ub;
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
    void SecondBoundaryConditions();
    //void SolveSLAE();
    void GeneratePortrait();
   

protected:
    vector<int> ia;
    vector<int> ja;
    vector<double> al;
    vector<double> di;
    vector<double> x;
    vector<double> b;
    vector<vector<int>> iaja;
    double eps = 1e-13;
    int maxIter = 10000;

    double Gauss3_Gxy(int i, int j, double b1, double b2, double b3, double b4, double b5, double b6, double a0, double a1, double a2);
    double Phi(double e, double n, int i, int j, double b1, double b2, double b3, double b4, double b5, double b6, double a0, double a1, double a2);
    double FUNteta(int number, double x, double y, double z);
    void AllocateMemory();
};

