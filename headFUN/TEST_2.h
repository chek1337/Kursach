#pragma once
#include "GRID.h"
double GridAndSLAE::FUN(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return x * y * z + 2 * x * y + 3 * y * z + 4 * x + 5 * y + 6;
		break;
	default:
		throw "������ � FUN";
		break;
	}
}

double GridAndSLAE::LAMBDA(int number)
{
	switch (number)
	{
	case 0:
		return 1;
		break;
	default:
		throw "������ � LAMBDA";
		break;
	}
}

double GridAndSLAE::GAMMA(int number)
{
	switch (number)
	{
	case 0:
		return 1;
		break;
	default:
		throw "������ � GAMMA";
		break;
	}
}



double GridAndSLAE::TETA(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return -(x * y + 3 * y);
		break;
	case 1:
		return (x * y + 3 * y);
			break;
	case 2:
		return -(y * z + 2 * y + 4);
			break;
	case 3:
		return (y * z + 2 * y + 4);
		break;
	case 4:
		return -(x*z + 2*x + 3*z + 5);
			break;
	case 5:
		return (x * z + 2 * x + 3 * z + 5);
			break;
	default:
		throw "������ � TETA";
		break;
	}
}

double GridAndSLAE::UBETA(int number, double x, double y, double z)
{
	switch (number)
	{
	default:
		throw "������ � UBETA";
		break;
	}
}