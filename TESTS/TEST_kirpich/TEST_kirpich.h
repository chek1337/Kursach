#pragma once
#include "GRID.h"
double GridAndSLAE::FUN(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return 0.4 * (5. + 0.2 * x + y + 30. * z + 0.5 * x * y + x * z + 10. * y * z + x * y * z);
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
		return 5.;
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
		return 0.4;
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
	case 0: // ƒл€ теста с 1им  Ё
		return -1. - 2.5 * y - 5. * z - 5. * y * z;
		break;
	case 1:
		return 5. + 2.5 * x + 50. * z + 5. * x * z;
		break;
	case 2:
		return 150. + 5. * x + 50. * y + 5. * x * y;
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
		return 9.1 + 11.25 * y + 50.5 * z + 30.5 * y * z;
		break;
	case 1:
		return 4.5 - 0.05 * x + 25. * z + 0.5 * x * z;
		break;
	case 2:
		return -10. - 0.3 * x - 4 * y;
	default:
		throw "ќшибка в UBETA";
		break;
	}
}