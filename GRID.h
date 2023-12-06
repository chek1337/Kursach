#pragma once
#include <iostream>
#include <vector>
#include <math.h>
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


class GridAndSLAE
{
public:
    // NumberOfNodes = NoN
    int NoN; // Общее кол-во узлов
    int NoN_xy; // Общее кол-во узлов в плоскости XY
    int NoN_z; // Общее кол-во узлов на оси Z
    int NoN_fe; // Общее кол-во конечных элементов на оси Z

    vector<FiniteElement> fe; // Массив содержащий информацию об о конечном элементе
    vector<Point> xy; // Массив с координатами сетки в Oxy
    vector<double> z; // Массив с координатами сетки вдоль оси Z
    void InputFromFile(FILE* inFE, FILE* inXY, FILE* inZ);

    void CalculateA_b();
    void SolveSLAE();

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
    void AllocateMemory();
};

