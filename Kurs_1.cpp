#include "GRID.h"

double FUN(int i, double x, double y, double z)
{
	switch (i)
	{
	case 0:
		return x * y * z;
		break;
	case 1:
		return sin(x * y * z);
		break;
	default:
		break;
	}
}


void GridAndSLAE::InputFromFile(FILE* inFE, FILE* inXY, FILE* inZ)
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
	//FormPortrait();
	GeneratePortrait();
	AllocateMemory();
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
	for (int curFE = 0; curFE < NoN_fe; curFE++)
	{
		vector <int> nodes_global;
		nodes_global.resize(8);
		nodes_global[0] = fe[curFE].node1 + fe[curFE].bottom * NoN_xy;
		nodes_global[1] = fe[curFE].node2 + fe[curFE].bottom * NoN_xy;
		nodes_global[2] = fe[curFE].node3 + fe[curFE].bottom * NoN_xy;
		nodes_global[3] = fe[curFE].node4 + fe[curFE].bottom * NoN_xy;

		double x1 = xy[nodes_global[0]].x;
		double y1 = xy[nodes_global[0]].y;

		double x2 = xy[nodes_global[1]].x;
		double y2 = xy[nodes_global[1]].y;

		double x3 = xy[nodes_global[2]].x;
		double y3 = xy[nodes_global[2]].y;

		double x4 = xy[nodes_global[3]].x;
		double y4 = xy[nodes_global[3]].y;

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
					M_local[i][j] = Mxy[mu_i][mu_j] * Mz[nu_i][nu_j]; // не забыть в дальнейшем про гамму
				}
			}

			for (int i = 0; i < 2; i++) // Формирование локальной жесткости для оси Z
			{
				for (int j = 0; j <= i; j++)
				{
					Gz[i][j] = Gz_0[i][j] / h_z;
				}
			}

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
					G_local[i][j] = Gxy[mu_i][mu_j] * Mz[nu_i][nu_j] + Mxy[mu_i][mu_j] * Gz[nu_i][nu_j]; // Не забыть потом про лямбду
				}
			}

			// Самое веселое начинается тут
			// Нужно закинуть элементы из лок. матриц в глоб матрицу
			// Так как я пока что забил на краевые условия, то
			// чтобы получить ij элемент матрицы А нужно Aij = Mij + Gij.
			// Номера узлов и соотв-их им базисных функций мы знаем из fe[i].nodex
			// То есть задача сводится к вытягиванию номера узла и сопаставлние ему
			// локальной базисной функции

			vector<double> b_local;
			b.resize(8);

			int numOfFun = 0; // нужен будет потом, когда будут разные функции на разных КЭ 
			vector<double> f_local;
			f_local.resize(8);

			f_local[0] = FUN(numOfFun, x1, y1, z1); // Вот здесь высока вероятность ошибки
			f_local[1] = FUN(numOfFun, x2, y2, z1); // так как я неявно считаю значения функции
			f_local[2] = FUN(numOfFun, x3, y3, z1); // в лок. узлах. Если бы нумерации узлов в КЭ поменялась
			f_local[3] = FUN(numOfFun, x4, y4, z1); // то может выскочить ошибка
			f_local[4] = FUN(numOfFun, x1, y1, z2);
			f_local[5] = FUN(numOfFun, x2, y2, z2);
			f_local[6] = FUN(numOfFun, x3, y3, z2);
			f_local[7] = FUN(numOfFun, x4, y4, z2);

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
				b[i] = sum;
			}

			double Aij, Aii;
			double lambda = 1, gamma = 1; //Не забыть их потом как кусочно лин. функциями сделать
			for (int i = 0; i < 8; i++)
			{
				Aii = gamma * M_local[i][i] + lambda * G_local[i][i];
				di[nodes_global[i]] += Aii;

				for (int j = i-1; j >=0; j--)
				{
					for (int k = ia[nodes_global[i]]; k < ia[nodes_global[i]+1]; k++)
					{ 
						if (ja[k] == nodes_global[j])
						{
							Aij = gamma * M_local[i][j] + lambda * G_local[i][j]; // Неправильно
							al[k] += Aij;
						}
					}
				}
			}

			nodes_global[0] = nodes_global[4];
			nodes_global[1] = nodes_global[5];
			nodes_global[2] = nodes_global[6];
			nodes_global[3] = nodes_global[7];
		}
		
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

int main()
{
	FILE* inFE, *inXY, *inZ;
	fopen_s(&inFE, "FE.txt", "r");
	fopen_s(&inXY, "XY.txt", "r");
	fopen_s(&inZ, "Z.txt", "r");
    GridAndSLAE grid;
	grid.InputFromFile(inFE, inXY, inZ);
	grid.CalculateA_b();
}
