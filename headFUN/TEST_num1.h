#pragma once
#include "C:\Users\Afina\source\repos\Kursach\Kurs_1\GRID.h"
double GridAndSLAE::FUN(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return x;
		break;
	default:
		throw "Ошибка в FUN";
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
		throw "Ошибка в LAMBDA";
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
		throw "Ошибка в GAMMA";
		break;
	}
}



double GridAndSLAE::TETA(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return 1;
		break;
	default:
		throw "Ошибка в TETA";
		break;
	}
}

double GridAndSLAE::UBETA(int number, double x, double y, double z)
{
	switch (number)
	{
	default:
		throw "Ошибка в UBETA";
		break;
	}
}