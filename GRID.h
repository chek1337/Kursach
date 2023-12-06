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
    int NoN; // ����� ���-�� �����
    int NoN_xy; // ����� ���-�� ����� � ��������� XY
    int NoN_z; // ����� ���-�� ����� �� ��� Z
    int NoN_fe; // ����� ���-�� �������� ��������� �� ��� Z

    vector<FiniteElement> fe; // ������ ���������� ���������� �� � �������� ��������
    vector<Point> xy; // ������ � ������������ ����� � Oxy
    vector<double> z; // ������ � ������������ ����� ����� ��� Z
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

