#pragma once
#include "GRID.h"
#include "math.h"
double GridAndSLAE::FUN(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return -(2 * y * cos(z) - x * x * y * cos(z)) + x * x * y * cos(z);
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