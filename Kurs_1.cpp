#include "GRID.h"

double FUN(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return 0.4*(5. + 0.2*x + y + 30.*z + 0.5*x*y + x*z + 10.*y*z + x*y*z);
		break;
	case 1:

		break;
	default:
		break;
	}
}

double LAMBDA(int number)
{
	switch (number)
	{
	case 0:
		return 5.;
		break;
	case 1:

		break;
	default:
		break;
	}
}

double GAMMA(int number)
{
	switch (number)
	{
	case 0:
		return 0.4;
		break;
	case 1:

		break;
	default:
		break;
	}
}



double GridAndSLAE::TETA(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return -1. - 2.5 * y - 5. * z - 5. * y * z;
		break;
	case 1:
		return 5. + 2.5 * x + 50. * z + 5. * x * z;
		break;
	case 2:
		return 150. + 5. * x + 50. * y + 5. * x * y;
		break;
	default:
		cout << "Ошибка в TETA";
		break;
	}
}

double GridAndSLAE::UBETA(int number, double x, double y, double z)
{
	switch (number)
	{
	case 0:
		return 9.1 + 11.25 * y + 50.5 * z + 30.5 * y * z;
		break;
	case 1:
		return 4.5 - 0.05 * x + 25. * z + 0.5 * x * z;
		break;
	case 2:
		return -10. - 0.3 * x - 4 * y;
		break;
	default:
		cout << "Ошибка в UBETA";
		break;
	}
}


void GridAndSLAE::InputFromFile(FILE* inFE, FILE* inXY, FILE* inZ, FILE* inFirstBC, FILE* inSecondBC, FILE* inThirdBC)
{
	// Запись информации об узлах конечных элементах
	fscanf_s(inFE, "%d", &NoN_fe);
	fe.resize(NoN_fe);
	// Мб тут сделать специально сортировку, чтобы номер узлов всегда шли бы на увелечение
	for (int i = 0; i < NoN_fe; i++)
	{
		fscanf_s(inFE, "%d", &fe[i].node1);
		fscanf_s(inFE, "%d", &fe[i].node2);
		fscanf_s(inFE, "%d", &fe[i].node3);
		fscanf_s(inFE, "%d", &fe[i].node4);
		fscanf_s(inFE, "%d", &fe[i].bottom);
		fscanf_s(inFE, "%d", &fe[i].top);
		fscanf_s(inFE, "%d", &fe[i].region);
	}

	// Запись информации об точках в Oxy
	fscanf_s(inXY, "%d", &NoN_xy);
	xy.resize(NoN_xy);
	for (int i = 0; i < NoN_xy; i++)
	{
		fscanf_s(inXY, "%lf", &xy[i].x);
		fscanf_s(inXY, "%lf", &xy[i].y);
	}

	// Запись информации об узлах вдоль оси Z
	fscanf_s(inZ, "%d", &NoN_z);
	z.resize(NoN_z);
	for (int i = 0; i < NoN_z; i++)
	{
		fscanf_s(inZ, "%lf", &z[i]);
	}
	NoN = NoN_xy * NoN_z;

	fscanf_s(inFirstBC, "%d", &Nof_FisrtBC);
	FirstBC.resize(Nof_FisrtBC);
	for (int i = 0; i < Nof_FisrtBC; i++)
	{
		fscanf_s(inFirstBC, "%d", &FirstBC[i].node);
		fscanf_s(inFirstBC, "%lf", &FirstBC[i].ug);
	}

	fscanf_s(inSecondBC, "%d", &Nof_SecondBC);
	SecondBC.resize(Nof_SecondBC);
	for (int i = 0; i < Nof_SecondBC; i++)
	{
		fscanf_s(inSecondBC, "%d", &SecondBC[i].node1);
		fscanf_s(inSecondBC, "%d", &SecondBC[i].node2);
		fscanf_s(inSecondBC, "%d", &SecondBC[i].node3);
		fscanf_s(inSecondBC, "%d", &SecondBC[i].node4);
		fscanf_s(inSecondBC, "%d", &SecondBC[i].num_teta);
	}

	fscanf_s(inThirdBC, "%d", &Nof_ThirdBC);
	ThirdBC.resize(Nof_ThirdBC);
	for (int i = 0; i < Nof_ThirdBC; i++)
	{
		fscanf_s(inThirdBC, "%d", &ThirdBC[i].node1);
		fscanf_s(inThirdBC, "%d", &ThirdBC[i].node2);
		fscanf_s(inThirdBC, "%d", &ThirdBC[i].node3);
		fscanf_s(inThirdBC, "%d", &ThirdBC[i].node4);
		fscanf_s(inThirdBC, "%d", &ThirdBC[i].num_ubeta);
		fscanf_s(inThirdBC, "%lf", &ThirdBC[i].beta);
	}
}

double sign(double x)
{
	if (x > 0)
		return 1;
	else if (x < 0)
		return -1;
	else
		return 0;
}

void GridAndSLAE::AllocateMemory()
{
	di.resize(NoN);
	iaja.resize(NoN);
	b.resize(NoN);
	x.resize(NoN);
	al.resize(ja.size());
}

void GridAndSLAE::CalculateA_b()
{
	GeneratePortrait();
	AllocateMemory();

	vector <int> nodes_global;
	nodes_global.resize(8);

	vector<double> b_local;
	b_local.resize(8);

	// нужен будет потом, когда будут разные функции на разных КЭ 
	vector<double> f_local;
	f_local.resize(8);
	//
	// Хранится только нижний треугольник
	//
	// Локальные матрицы массы
	//
	vector<vector<double>> M0 = { {4}, {2, 4}, {2, 1, 4}, {1, 2, 2, 4} };
	vector<vector<double>> M1 = { {2}, {2, 6}, {1, 1, 2}, {1, 2, 3, 6} };
	vector<vector<double>>  M2 = { {2}, {1, 2}, {2, 1, 6}, {1, 2, 3, 6} };
	vector<vector<double>> Mxy;
	Mxy.resize(4);
	for (int i = 0; i < 4; i++)
		Mxy[i].resize(i + 1);

	vector<vector<double>>  Mz_0 = { {2}, {1, 2} };
	vector<vector<double>>  Mz;
	Mz.resize(2);
	for (int i = 0; i < 2; i++)
		Mz[i].resize(i + 1);
	

	vector<vector<double>> M_local;
	M_local.resize(8); 
	for (int i = 0; i < 8; i++)
		M_local[i].resize(i + 1);
	

	//
	// Локальые матрицы жесткости
	//
	vector<vector<double>> Gxy;
	Gxy.resize(4);
	for (int i = 0; i < 4; i++)
		Gxy[i].resize(i + 1);
	

	vector<vector<double>> Gz_0 = { {1}, {-1, 1} };
	vector<vector<double>> Gz;
	Gz.resize(2);
	for (int i = 0; i < 2; i++)
		Gz[i].resize(i + 1);
	
	vector<vector<double>> G_local;
	G_local.resize(8);
	for (int i = 0; i < 8; i++)
		G_local[i].resize(i + 1);

	//
	// Сборка локальных матриц
	//
	for (int curFE = 0; curFE < NoN_fe; )
	{
		
		nodes_global[0] = fe[curFE].node1 + fe[curFE].bottom * NoN_xy;
		nodes_global[1] = fe[curFE].node2 + fe[curFE].bottom * NoN_xy;
		nodes_global[2] = fe[curFE].node3 + fe[curFE].bottom * NoN_xy;
		nodes_global[3] = fe[curFE].node4 + fe[curFE].bottom * NoN_xy;

		double region_cur = fe[curFE].region;

		double x1 = xy[fe[curFE].node1].x;
		double y1 = xy[fe[curFE].node1].y;

		double x2 = xy[fe[curFE].node2].x;
		double y2 = xy[fe[curFE].node2].y;

		double x3 = xy[fe[curFE].node3].x;
		double y3 = xy[fe[curFE].node3].y;

		double x4 = xy[fe[curFE].node4].x;
		double y4 = xy[fe[curFE].node4].y;

		double a0 = ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1));
		double a1 = ((x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3));
		double a2 = ((y3 - y1) * (x4 - x2) - (x3 - x1) * (y4 - y2));

		for (int i = 0; i < 4; i++) // Формирование локальной матрицы массы для Oxy
			for (int j = 0; j <= i; j++)
				Mxy[i][j] = sign(a0) * (a0 / 36. * M0[i][j] + a1 / 72. * M1[i][j] + a2 / 72. * M2[i][j]);

		double b1 = x3 - x1;
		double b2 = x2 - x1;
		double b5 = x1 - x2 - x3 + x4;

		double b3 = y3 - y1;
		double b4 = y2 - y1;
		double b6 = y1 - y2 - y3 + y4;

		for (int i = 0; i < 4; i++)
			for (int j = 0; j <= i; j++) // Формирование локальной жесткости для Oxy
				Gxy[i][j] = Gauss3_Gxy(i, j, b1, b2, b3, b4, b5, b6, a0, a1, a2);

		// Короче, поменял обход области (последовательнсти КЭ). До этого было так, я обходил сначала все конечные элементы,
		// у которых высота 0-1, потом 1-2, потом 2-3, и тд
		// Но я тут заметил, что по сути КЭ, которые "стоят" в одном столбце имеют одинаковое основание
		// Поэтому достотаточно один раз посчитать Mxy, Gxy, а Mz, Gz, пересчитывать на каждом уровне

		for (int lvl = 0; lvl < NoN_z - 1; lvl++, curFE++) // Здесь как раз за это и отвечает этот цикл
		{

			nodes_global[4] = fe[curFE].node1 + fe[curFE].top * NoN_xy;
			nodes_global[5] = fe[curFE].node2 + fe[curFE].top * NoN_xy;
			nodes_global[6] = fe[curFE].node3 + fe[curFE].top * NoN_xy;
			nodes_global[7] = fe[curFE].node4 + fe[curFE].top * NoN_xy;

			double z1 = z[fe[curFE].bottom];
			double z2 = z[fe[curFE].top];

			double h_z = z2 - z1;

			for (int i = 0; i < 2; i++) // Формирование локальной матрицы массы для оси Z
				for (int j = 0; j <= i; j++)
					Mz[i][j] = h_z / 6. * Mz_0[i][j];

			for (int i = 0; i < 8; i++) // Сборка локальной матрицы массы
			{
				for (int j = 0; j <= i; j++)
				{
					int mu_i = i % 4;
					int mu_j = j % 4;
					if (mu_i < mu_j)
						swap(mu_i, mu_j);
					int nu_i = (int)(i / 4.);
					int nu_j = (int)(j / 4.);
					/*if (nu_i < nu_j)
						swap(nu_i, nu_j);*/
					M_local[i][j] = Mxy[mu_i][mu_j] * Mz[nu_i][nu_j];
				}
			}

			for (int i = 0; i < 2; i++) // Формирование локальной жесткости для оси Z
				for (int j = 0; j <= i; j++)
					Gz[i][j] = Gz_0[i][j] / h_z;

			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j <= i; j++) // Сборка локальной матрицы жесткости
				{
					int mu_i = i % 4;
					int mu_j = j % 4;
					if (mu_i < mu_j)
						swap(mu_i, mu_j);
					int nu_i = (int)(i / 4.);
					int nu_j = (int)(j / 4.);
					G_local[i][j] = Gxy[mu_i][mu_j] * Mz[nu_i][nu_j] + Mxy[mu_i][mu_j] * Gz[nu_i][nu_j];
				}
			}

			// Самое веселое начинается тут
			// Нужно закинуть элементы из лок. матриц в глоб матрицу
			// Так как я пока что забил на краевые условия, то
			// чтобы получить ij элемент матрицы А нужно Aij = Mij + Gij.
			// Номера узлов и соотв-их им базисных функций мы знаем из fe[i].nodex
			// То есть задача сводится к вытягиванию номера узла и сопаставлние ему
			// локальной базисной функции

			f_local[0] = FUN(region_cur, x1, y1, z1); // Вот здесь высока вероятность ошибки
			f_local[1] = FUN(region_cur, x2, y2, z1); // так как я неявно считаю значения функции
			f_local[2] = FUN(region_cur, x3, y3, z1); // в лок. узлах. Если бы нумерации узлов в КЭ поменялась
			f_local[3] = FUN(region_cur, x4, y4, z1); // то может выскочить ошибка
			f_local[4] = FUN(region_cur, x1, y1, z2);
			f_local[5] = FUN(region_cur, x2, y2, z2);
			f_local[6] = FUN(region_cur, x3, y3, z2);
			f_local[7] = FUN(region_cur, x4, y4, z2);

			for (int i = 0; i < 8; i++) //Нужно потом более оптимально  умножение сделать
			{
				double sum = 0;
				for (int j = 0; j < 8; j++)
				{
					if (i > j)
						sum += f_local[j] * M_local[i][j];
					else
						sum += f_local[j] * M_local[j][i];
				}
				b_local[i] = sum;
			}

			double Aij, Aii;
			double lambda = LAMBDA(region_cur);
			double gamma = GAMMA(region_cur);
			for (int i = 0; i < 8; i++)
			{
				Aii = gamma * M_local[i][i] + lambda * G_local[i][i];
				di[nodes_global[i]] += Aii;

				for (int j = i-1; j >= 0; j--)
				{
					for (int k = ia[nodes_global[i]]; k < ia[nodes_global[i]+1]; k++)
					{ 
						if (ja[k] == nodes_global[j])
						{
							Aij = gamma * M_local[i][j] + lambda * G_local[i][j];
							al[k] += Aij;
						}
					}
				}
			}

			for (int i = 0; i < 8; i++)
				b[nodes_global[i]] += b_local[i];

			nodes_global[0] = nodes_global[4];
			nodes_global[1] = nodes_global[5];
			nodes_global[2] = nodes_global[6];
			nodes_global[3] = nodes_global[7];
		}
		
	}
}

void GridAndSLAE::SecondBoundaryConditions()
{
	vector<double> bS2_local(4, 0);
	vector<int> nodes_global(4, 0);
	vector<double> teta_local(4, 0);

	//
	// Матрицы для боковой грани
	//

	vector<vector<double>>  MXorY_0{ {2}, {1, 2} };
	vector<vector<double>>  MXorY;
	MXorY.resize(2);
	for (int i = 0; i < 2; i++)
		MXorY[i].resize(i + 1);

	vector<vector<double>>  Mz_0 = { {2}, {1, 2} };
	vector<vector<double>>  Mz;
	Mz.resize(2);
	for (int i = 0; i < 2; i++)
		Mz[i].resize(i + 1);

	vector<vector<double>>  M_XYorXZorYZ; 
	M_XYorXZorYZ.resize(4);
	for (int i = 0; i < 4; i++)
		M_XYorXZorYZ[i].resize(i + 1);
	// 
	//******************************************
	// 

	// 
	// Матрицы для основания
	// 
	vector<vector<double>> M0 = { {4}, {2, 4}, {2, 1, 4}, {1, 2, 2, 4} };
	vector<vector<double>> M1 = { {2}, {2, 6}, {1, 1, 2}, {1, 2, 3, 6} };
	vector<vector<double>>  M2 = { {2}, {1, 2}, {2, 1, 6}, {1, 2, 3, 6} };
	// 
	//******************************************
	//

	for (int curSecondBC = 0; curSecondBC < Nof_SecondBC; curSecondBC++)
	{
		
		nodes_global[0] = SecondBC[curSecondBC].node1;
		nodes_global[1] = SecondBC[curSecondBC].node2;
		nodes_global[2] = SecondBC[curSecondBC].node3;
		nodes_global[3] = SecondBC[curSecondBC].node4;
		int num_teta_local = SecondBC[curSecondBC].num_teta;

		// Определить что это - боковая грань или основание? 
		if ((nodes_global[0] % NoN_xy) == (nodes_global[2] % NoN_xy) && (nodes_global[1] % NoN_xy) == (nodes_global[3] % NoN_xy))
		{
			// Значит боковая грань
			// Будем обозначать матрицу массы для x/y как M_XorY (смотря паралельно какой оси находится грань)
			// Будем искать матрицу массы для краевого условия как M_XorY * Mz

			double x1 = xy[nodes_global[0] % NoN_xy].x;
			double x2 = xy[nodes_global[1] % NoN_xy].x;
			double y1 = xy[nodes_global[0] % NoN_xy].y;
			double y2 = xy[nodes_global[1] % NoN_xy].y;
			double h_XorY = 0;

			div_t result = div(nodes_global[0], NoN_xy);
			double z1 = z[result.quot];
			double z2 = z[result.quot + 1];
			double h_z = z2 - z1;
			
			teta_local[0] = TETA(num_teta_local, x1, y1, z1);
			teta_local[1] = TETA(num_teta_local, x2, y2, z1);
			teta_local[2] = TETA(num_teta_local, x1, y1, z2);
			teta_local[3] = TETA(num_teta_local, x2, y2, z2);

			if ((nodes_global[0] + 1) == nodes_global[1]) // значит боковая грань поралельна Ox. Я предполагаю что нумерация у узлов корректная
				h_XorY = x2 - x1;
			else if ((nodes_global[0] % NoN_xy) == (nodes_global[2] % NoN_xy)) // значит боковая грань поралельна Oy.
				h_XorY = y2 - y1;
			else 
				cout << "1 Словил какую-ту херь в 2ом краевом\n";

			for (int i = 0; i < 2; i++)
				for (int j = 0; j <= i; j++)
					MXorY[i][j] = h_XorY / 6. * MXorY_0[i][j];

			for (int i = 0; i < 2; i++)
				for (int j = 0; j <= i; j++)
					Mz[i][j] = h_z / 6. * Mz_0[i][j];

			// стр 234 кирпича
			// Вообще можно сделать сразу матрицу MXZorYZ размером 4х4 (стр 233 5.20), мб потом переделаю
			// если кто-то придумает как это красиво можно в цикле сделать - буду рад
			// а пока тупо скатаю формулы (при учете того, что у меня нижний треугольник, а в кирпиче верхний)
			// но лучше перепроверить все равно
			M_XYorXZorYZ[0][0] = MXorY[0][0] * Mz[0][0];
			M_XYorXZorYZ[1][0] = MXorY[1][0] * Mz[0][0];
			M_XYorXZorYZ[2][0] = MXorY[0][0] * Mz[1][0];
			M_XYorXZorYZ[3][0] = MXorY[1][0] * Mz[1][0];

			M_XYorXZorYZ[1][1] = MXorY[1][1] * Mz[0][0];
			M_XYorXZorYZ[2][1] = MXorY[1][0] * Mz[1][0];
			M_XYorXZorYZ[3][1] = MXorY[1][1] * Mz[1][0];

			M_XYorXZorYZ[2][2] = MXorY[0][0] * Mz[1][1];
			M_XYorXZorYZ[3][2] = MXorY[1][0] * Mz[1][1];

			M_XYorXZorYZ[3][3] = MXorY[1][1] * Mz[1][1];

		}
		else if((nodes_global[0] + 1) == nodes_global[1] && (nodes_global[2] + 1) == nodes_global[3])
		{
			// Иначе основание

			double x1 = xy[nodes_global[0] % NoN_xy].x;
			double y1 = xy[nodes_global[0] % NoN_xy].y;

			double x2 = xy[nodes_global[1] % NoN_xy].x;
			double y2 = xy[nodes_global[1] % NoN_xy].y;

			double x3 = xy[nodes_global[2] % NoN_xy].x;
			double y3 = xy[nodes_global[2] % NoN_xy].y;

			double x4 = xy[nodes_global[3] % NoN_xy].x;
			double y4 = xy[nodes_global[3] % NoN_xy].y;

			div_t result = div(nodes_global[0], NoN_xy);
			double zlvl = z[result.quot];

			double a0 = ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1));
			double a1 = ((x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3));
			double a2 = ((y3 - y1) * (x4 - x2) - (x3 - x1) * (y4 - y2));

			teta_local[0] = TETA(num_teta_local, x1, y1, zlvl);
			teta_local[1] = TETA(num_teta_local, x2, y2, zlvl);
			teta_local[2] = TETA(num_teta_local, x3, y3, zlvl);
			teta_local[3] = TETA(num_teta_local, x4, y4, zlvl);

			for (int i = 0; i < 4; i++) // Формирование локальной матрицы массы для Oxy
				for (int j = 0; j <= i; j++)
					M_XYorXZorYZ[i][j] = sign(a0) * (a0 / 36. * M0[i][j] + a1 / 72. * M1[i][j] + a2 / 72. * M2[i][j]);
		}
		else
			"2 Словил какую-ту херь в 2ом краевом\n";

		for (int i = 0; i < 4; i++) //Нужно потом более оптимально  умножение сделать
		{
			double sum = 0;
			for (int j = 0; j < 4; j++)
			{
				if (i > j)
					sum += teta_local[j] * M_XYorXZorYZ[i][j];
				else
					sum += teta_local[j] * M_XYorXZorYZ[j][i];
			}
			bS2_local[i] = sum;
		}
		for (int i = 0; i < 4; i++)
			b[nodes_global[i]] += bS2_local[i];
	}
}



void GridAndSLAE::ThirdBoundaryConditions() // ctrl+c -> ctrl+v из SecondBoundaryConditions. 
{											// Так что если и есть ошибки, то они могли всплыть из-за невдумчивого копирования 
	vector<double> bS3_local(4, 0);
	vector<int> nodes_global(4, 0);
	vector<double> ubeta_local(4, 0);

	//
	// Матрицы для боковой грани
	//
	vector<vector<double>>  MXorY_0{ {2}, {1, 2} };
	vector<vector<double>>  MXorY;
	MXorY.resize(2);
	for (int i = 0; i < 2; i++)
		MXorY[i].resize(i + 1);

	vector<vector<double>>  Mz_0 = { {2}, {1, 2} };
	vector<vector<double>>  Mz;
	Mz.resize(2);
	for (int i = 0; i < 2; i++)
		Mz[i].resize(i + 1);

	vector<vector<double>>  M_XYorXZorYZ;
	M_XYorXZorYZ.resize(4);
	for (int i = 0; i < 4; i++)
		M_XYorXZorYZ[i].resize(i + 1);

	// 
	// Матрицы для основания
	// 
	vector<vector<double>> M0 = { {4}, {2, 4}, {2, 1, 4}, {1, 2, 2, 4} };
	vector<vector<double>> M1 = { {2}, {2, 6}, {1, 1, 2}, {1, 2, 3, 6} };
	vector<vector<double>>  M2 = { {2}, {1, 2}, {2, 1, 6}, {1, 2, 3, 6} };
	// 
	//******************************************
	//
	for (int curThirdBC = 0; curThirdBC < Nof_ThirdBC; curThirdBC++)
	{

		nodes_global[0] = ThirdBC[curThirdBC].node1;
		nodes_global[1] = ThirdBC[curThirdBC].node2;
		nodes_global[2] = ThirdBC[curThirdBC].node3;
		nodes_global[3] = ThirdBC[curThirdBC].node4;

		int num_ubeta_local = ThirdBC[curThirdBC].num_ubeta;
		// Определить что это - боковая грань или основание? 
		if ((nodes_global[0] % NoN_xy) == (nodes_global[2] % NoN_xy) && (nodes_global[1] % NoN_xy) == (nodes_global[3] % NoN_xy))
		{
			// Значит боковая грань
			// Будем обозначать матрицу массы для x/y как M_XorY (смотря паралельно какой оси находится грань)
			// Будем искать матрицу массы для краевого условия как M_XorY * Mz

			double x1 = xy[nodes_global[0] % NoN_xy].x;
			double x2 = xy[nodes_global[1] % NoN_xy].x;
			double y1 = xy[nodes_global[0] % NoN_xy].y;
			double y2 = xy[nodes_global[1] % NoN_xy].y;
			double h_XorY = 0;

			div_t result = div(nodes_global[0], NoN_xy);
			double z1 = z[result.quot];
			double z2 = z[result.quot + 1];
			double h_z = z2 - z1;

			ubeta_local[0] = TETA(num_ubeta_local, x1, y1, z1);
			ubeta_local[1] = TETA(num_ubeta_local, x2, y2, z1);
			ubeta_local[2] = TETA(num_ubeta_local, x1, y1, z2);
			ubeta_local[3] = TETA(num_ubeta_local, x2, y2, z2);

			if ((nodes_global[0] + 1) == nodes_global[1]) // значит боковая грань поралельна Ox. Я предполагаю что нумерация у узлов корректная
				h_XorY = x2 - x1;
			else if ((nodes_global[0] % NoN_xy) == (nodes_global[2] % NoN_xy)) // значит боковая грань поралельна Oy.
				h_XorY = y2 - y1;
			else
				cout << "1 Словил какую-ту херь в 3ем краевом\n";

			for (int i = 0; i < 2; i++)
				for (int j = 0; j <= i; j++)
					MXorY[i][j] = h_XorY / 6. * MXorY_0[i][j];

			for (int i = 0; i < 2; i++)
				for (int j = 0; j <= i; j++)
					Mz[i][j] = h_z / 6. * Mz_0[i][j];

			// стр 234 кирпича
			// Вообще можно сделать сразу матрицу MXZorYZ размером 4х4 (стр 233 5.20), мб потом переделаю
			// если кто-то придумает как это красиво можно в цикле сделать - буду рад
			// а пока тупо скатаю формулы (при учете того, что у меня нижний треугольник, а в кирпиче верхний)
			// но лучше перепроверить все равно
			M_XYorXZorYZ[0][0] = MXorY[0][0] * Mz[0][0];
			M_XYorXZorYZ[1][0] = MXorY[1][0] * Mz[0][0];
			M_XYorXZorYZ[2][0] = MXorY[0][0] * Mz[1][0];
			M_XYorXZorYZ[3][0] = MXorY[1][0] * Mz[1][0];

			M_XYorXZorYZ[1][1] = MXorY[1][1] * Mz[0][0];
			M_XYorXZorYZ[2][1] = MXorY[1][0] * Mz[1][0];
			M_XYorXZorYZ[3][1] = MXorY[1][1] * Mz[1][0];

			M_XYorXZorYZ[2][2] = MXorY[0][0] * Mz[1][1];
			M_XYorXZorYZ[3][2] = MXorY[1][0] * Mz[1][1];

			M_XYorXZorYZ[3][3] = MXorY[1][1] * Mz[1][1];
		}
		else if ((nodes_global[0] + 1) == nodes_global[1] && (nodes_global[2] + 1) == nodes_global[3])
		{
			// Иначе основание

			double x1 = xy[nodes_global[0] % NoN_xy].x;
			double y1 = xy[nodes_global[0] % NoN_xy].y;

			double x2 = xy[nodes_global[1] % NoN_xy].x;
			double y2 = xy[nodes_global[1] % NoN_xy].y;

			double x3 = xy[nodes_global[2] % NoN_xy].x;
			double y3 = xy[nodes_global[2] % NoN_xy].y;

			double x4 = xy[nodes_global[3] % NoN_xy].x;
			double y4 = xy[nodes_global[3] % NoN_xy].y;

			div_t result = div(nodes_global[0], NoN_xy);
			double zlvl = z[result.quot];

			double a0 = ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1));
			double a1 = ((x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3));
			double a2 = ((y3 - y1) * (x4 - x2) - (x3 - x1) * (y4 - y2));

			ubeta_local[0] = TETA(num_ubeta_local, x1, y1, zlvl);
			ubeta_local[1] = TETA(num_ubeta_local, x2, y2, zlvl);
			ubeta_local[2] = TETA(num_ubeta_local, x3, y3, zlvl);
			ubeta_local[3] = TETA(num_ubeta_local, x4, y4, zlvl);

			for (int i = 0; i < 4; i++) // Формирование локальной матрицы массы для Oxy
				for (int j = 0; j <= i; j++)
					M_XYorXZorYZ[i][j] = sign(a0) * (a0 / 36. * M0[i][j] + a1 / 72. * M1[i][j] + a2 / 72. * M2[i][j]);
		}
		else
			"2 Словил какую-ту херь в 3ем краевом\n";

		// Добавка в глоб матрицу
		double beta = ThirdBC[curThirdBC].beta;
		for (int i = 0; i < 4; i++)
		{
			di[nodes_global[i]] += beta * M_XYorXZorYZ[i][i];
			for (int j = i - 1; j >= 0; j--)
				for (int k = ia[nodes_global[i]]; k < ia[nodes_global[i] + 1]; k++)
					if (ja[k] == nodes_global[j])
						al[k] += beta * M_XYorXZorYZ[i][j];
		}

		for (int i = 0; i < 4; i++) //Нужно потом более оптимально  умножение сделать
		{
			double sum = 0;
			for (int j = 0; j < 4; j++)
			{
				if (i > j)
					sum += ubeta_local[j] * M_XYorXZorYZ[i][j];
				else
					sum += ubeta_local[j] * M_XYorXZorYZ[j][i];
			}
			bS3_local[i] = sum;
		}
		for (int i = 0; i < 4; i++)
			b[nodes_global[i]] += bS3_local[i];
	}
}

void GridAndSLAE::GeneratePortrait()
{
	iaja.resize(NoN);
	vector<int> nodes_local(8, 0);
	for (int curFE = 0; curFE < NoN_fe; curFE++)
	{
		nodes_local[0] = fe[curFE].node1 + fe[curFE].bottom * NoN_xy;
		nodes_local[1] = fe[curFE].node2 + fe[curFE].bottom * NoN_xy;
		nodes_local[2] = fe[curFE].node3 + fe[curFE].bottom * NoN_xy;
		nodes_local[3] = fe[curFE].node4 + fe[curFE].bottom * NoN_xy;
		nodes_local[4] = fe[curFE].node1 + fe[curFE].top * NoN_xy;
		nodes_local[5] = fe[curFE].node2 + fe[curFE].top * NoN_xy;
		nodes_local[6] = fe[curFE].node3 + fe[curFE].top * NoN_xy;
		nodes_local[7] = fe[curFE].node4 + fe[curFE].top * NoN_xy;
		for (int i = 0; i < 8; i++)
		{
			int node_i = nodes_local[i];
			for (int j = 0; j < 8; j++)
			{
				int node_j = nodes_local[j];
				if (node_i > node_j)
				{
					int flag = 0;
					for (int k = 0; k < iaja[node_i].size() and flag != 1; k++)
					{
						if (iaja[node_i][k] == node_j) 
						{
							flag = 1;
							continue;
						}
					}
					if (!flag)
						iaja[node_i].push_back(node_j);
				}
			}
		}
	}
	// Сортировка
	for (int i = 0; i < NoN; i++)
		sort(iaja[i].begin(), iaja[i].end());

	ia.resize(NoN + 1);
	ia[0] = 0;
	for (int i = 0; i < NoN; i++)
	{
		int size = iaja[i].size();
		ia[i + 1] = ia[i] + size;
		for (int j = 0; j < size; j++)
			ja.push_back(iaja[i][j]);
	}
}

double GridAndSLAE::Phi(double e, double n, int i, int j, double b1, double b2, double b3, double b4, double b5, double b6, double a0, double a1, double a2) // Ммммммм список аргументов
{
	// Формула интеграла
	// SS (phi_i_1 * phi_j_1 + phi_i_2 * phi_j_2)/J dedn

	double J = 0;
	double value = 0;
	double phi_i_1 = 0, phi_j_1 = 0, phi_i_2 = 0, phi_j_2 = 0;

	switch (i) // Не стал делать для всяких "b6*e + b3" отдельные переменные, чтобы легче можно было формулам проверять 
	{
	case 0:
		// phi = (1 - e)*(1 - n)
		// d(phi)/de = -1 + n
		// d(phi)/dn = -1 + e
		phi_i_1 = (n - 1.) * (b6 * e + b3) - (e - 1.) * (b6 * n + b4);
		phi_i_2 = (e - 1.) * (b5 * n + b2) - (n - 1.) * (b5 * e + b1);
		break;
	case 1:
		// phi = e*(1 - n)
		// d(phi)/de = 1 - n
		// d(phi)/dn = -e
		phi_i_1 = (1. - n) * (b6 * e + b3) - (-e) * (b6 * n + b4);
		phi_i_2 = (-e) * (b5 * n + b2) - (1. - n) * (b5 * e + b1);
		break;
	case 2:
		// phi = (1 - e)*n
		// d(phi)/de = -n
		// d(phi)/dn = 1 - e
		phi_i_1 = (-n) * (b6 * e + b3) - (1. - e) * (b6 * n + b4);
		phi_i_2 = (1. - e) * (b5 * n + b2) - (-n) * (b5 * e + b1);
		break;
	case 3:
		// phi = e*n
		// d(phi)/de = n
		// d(phi)/dn = e
		phi_i_1 = n * (b6 * e + b3) - e * (b6 * n + b4);
		phi_i_2 = e * (b5 * n + b2) - n * (b5 * e + b1);
		break;
	}

	switch (j)
	{
	case 0:
		// phi = (1 - e)*(1 - n)
		// d(phi)/de = -1 + n
		// d(phi)/dn = -1 + e
		phi_j_1 = (n - 1.) * (b6 * e + b3) - (e - 1.) * (b6 * n + b4);
		phi_j_2 = (e - 1.) * (b5 * n + b2) - (n - 1.) * (b5 * e + b1);
		break;
	case 1:
		// phi = e*(1 - n)
		// d(phi)/de = 1 - n
		// d(phi)/dn = -e
		phi_j_1 = (1. - n) * (b6 * e + b3) - (-e) * (b6 * n + b4);
		phi_j_2 = (-e) * (b5 * n + b2) - (1. - n) * (b5 * e + b1);
		break;
	case 2:
		// phi = (1 - e)*n
		// d(phi)/de = -n
		// d(phi)/dn = 1 - e
		phi_j_1 = (-n) * (b6 * e + b3) - (1. - e) * (b6 * n + b4);
		phi_j_2 = (1. - e) * (b5 * n + b2) - (-n) * (b5 * e + b1);
		break;
	case 3:
		// phi = e*n
		// d(phi)/de = n
		// d(phi)/dn = e
		phi_j_1 = n * (b6 * e + b3) - e * (b6 * n + b4);
		phi_j_2 = e * (b5 * n + b2) - n * (b5 * e + b1);
		break;
	}
	J = a0 + a1 * e + a2 * n;
	value = (phi_i_1 * phi_j_1 + phi_i_2 * phi_j_2)/J;

	return value;
}

double GridAndSLAE::Gauss3_Gxy(int i, int j, double b1, double b2, double b3, double b4, double b5, double b6, double a0, double a1, double a2)
{

	double sum = 0;
	double tauK = 0, tauL = 0;
	double tK = 0, tL = 0;
	double q = sqrt(3. / 5.);

	sum += (8. / 9.) * (8. / 9.) * Phi((1. + 0) / 2., (1. + 0) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);
	sum += (8. / 9.) * (5. / 9.) * Phi((1. + 0) / 2., (1. + q) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);
	sum += (8. / 9.) * (5. / 9.) * Phi((1. + 0) / 2., (1. - q) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);

	sum += (5. / 9.) * (8. / 9.) * Phi((1. + q) / 2., (1. + 0) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);
	sum += (5. / 9.) * (5. / 9.) * Phi((1. + q) / 2., (1. + q) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);
	sum += (5. / 9.) * (5. / 9.) * Phi((1. + q) / 2., (1. - q) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);

	sum += (5. / 9.) * (8. / 9.) * Phi((1. - q) / 2., (1. + 0) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);
	sum += (5. / 9.) * (5. / 9.) * Phi((1. - q) / 2., (1. + q) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);
	sum += (5. / 9.) * (5. / 9.) * Phi((1. - q) / 2., (1. - q) / 2., /**/ i, j, /**/ b1, b2, b3, b4, b5, b6, a0, a1, a2);

	sum /= 4.;
	return sum;
}

void GridAndSLAE::OutputDense()
{
	int flagfound = 0;
	for (int i = 0; i < NoN; i++)
	{
		int k = ia[i + 1] - ia[i];
		if (k == 0)
		{
			for (int j = 0; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}
		else
		{
			int lastj = 0;
			for (int j = ia[i]; j < ia[i + 1]; j++) 
			{
				for (int p = lastj; p < ja[j]; p++) 
				{
					printf(REALOUTD, 0.0);
				}
				printf(REALOUTD, 216*al[j]); // ТУТ ДОБАВИЛ 216
				lastj = ja[j] + 1;
			}
			for (int j = lastj; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}

		printf(REALOUTD, 216 * di[i]); // ТУТ ДОБАВИЛ 216

		for (int j = i + 1; j < NoN; j++)
		{
			k = ia[j + 1] - ia[j];
			if (k == 0) {
				printf(REALOUTD, 0.0);
			}
			else
			{
				flagfound = 0;
				for (k = ia[j]; k < ia[j + 1]; k++)
				{

					if (ja[k] == i)
					{
						printf(REALOUTD, 216 * al[k]); // ТУТ ДОБАВИЛ 216
						flagfound = 1;
						break;
					}
				}
				if (flagfound == 0)
					printf(REALOUTD, 0.0);
			}
		}
		printf("\n");
	}
	printf("\n");

	for (int i = 0; i < NoN; i++)
	{
		printf("%.15lf\n", 216 * b[i]); // ТУТ ДОБАВИЛ 216
	}
}


