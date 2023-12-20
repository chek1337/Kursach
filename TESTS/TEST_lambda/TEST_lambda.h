#pragma once
#include "GRID.h"
double GridAndSLAE::FUN(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return 8 * x * (y * z + z + y + 1);
		break;
	case 1:
		return 3 * x * z * (y + 1);
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
		return 2./3.;
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
		return 4.;
		break;
	case 1:
		return 1.;
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
		return -2 * x * (z + 1);;
		break;
	case 1:
		return -2*x*z;
		break;
	case 2:
		return 2 * (y * z + z + y + 1);
		break;
	case 3:
		return 2*z*(1+y);
		break;
	case 4:
		return 2*x*(1+y);
	case 5:
		return -2*x*(1+y);
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
		return 8 * x * (z + 1);
		break;
	case 1:
		return 11.5 * x * z;
		break;
	default:
		throw "ќшибка в UBETA";
		break;
	}
}