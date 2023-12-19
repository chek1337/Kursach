#pragma once
#include "GRID.h"
double GridAndSLAE::FUN(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return 2 * (x * y * z + 2 * x * y + 3 * y * z + 4);
		break;
	case 1:
		return 3 * 0.1 * (x * y * z + 2 * x * y + 3 * y * z);
		break;

	default:
		throw "ќшибка в FUN";
		break;
	}
}

double GridAndSLAE::LAMBDA(int number)
{
	switch (number)
	{
	case 0:
		return 1.;
		break;
	case 1:
		return 10.;
		break;
	default:
		throw "ќшибка в LAMBDA";
		break;
	}
}

double GridAndSLAE::GAMMA(int number)
{
	switch (number)
	{
	case 0:
		return 2.;
		break;
	case 1:
		return 3.;
		break;
	default:
		throw "ќшибка в GAMMA";
		break;
	}
}



double GridAndSLAE::TETA(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return -(y*z+2*y);
		break;
	case 1:
		return x*z+2*x+3*z;
		break;
	case 2:
		return x * z + 2 * x + 3 * z;
	case 3:
		return x * y + 3 * y;
		break;
	case 4:
		return x * y + 3 * y;
		break;
	default:
		throw "ќшибка в TETA";
		break;
	}
}

double GridAndSLAE::UBETA(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:// ƒл€ теста с 1им  Ё
		return 0.1 * (x * y * z + 2 * x * y + 3 * y * z) + (y * z + 2 * y) / 2;
		break;
	case 1:
		return x * y * z + 2 * x * y + 3 * y * z + 4 - (x * z + 2 * x + 3 * z) / 2;
		break;
	case 2:
		return 0.1 * (x * y * z + 2 * x * y + 3 * y * z) - (x * z + 2 * x + 3 * z) / 2;
	case 3:
		return  x * y * z + 2 * x * y + 3 * y * z + 4 - (x * y + 3 * y) / 2;
		break;
	case 4:
		return  0.1*(x * y * z + 2 * x * y + 3 * y * z + 4) - (x * y + 3 * y) / 2;
		break;
	default:
		throw "ќшибка в UBETA";
		break;
	}
}