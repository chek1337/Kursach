#pragma once
#include "GRID.h"
#include "math.h"
double GridAndSLAE::FUN(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return x*x*x + y*y*y + z*z*z - (6*x + 6*y + 6*z);
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