#include "BandMatrix.h"
#include <vector>

using namespace std;

void BandMatrix::resize(int size)
{
	diag.resize(size);
	top.resize(size);
	bot.resize(size);
}

void BandMatrix::toZero()
{
	fill(diag.begin(), diag.end(), 0);
	fill(top.begin(), top.end(), 0);
	fill(bot.begin(), bot.end(), 0);
}

void BandMatrix::matrixToLU()
{
	for (int i = 1; i < diag.size(); i++)
	{
		bot[i] = bot[i] / diag[i - 1];
		diag[i] = diag[i] - bot[i] * top[i - 1];
	}
}

void BandMatrix::solveLUx(vector<double>& b, vector<double>& x)
{
	solveLy(b, x);
	solveUx(x, x);
}

void BandMatrix::solveUx(vector<double>& b, vector<double>& x)
{
	x[x.size() - 1] = b[b.size() - 1];
	for (int i = x.size() - 2; i >= 0; i--)
	{
		x[i] = (b[i] - top[i] * x[i + 1]) / diag[i];
	}
}

void BandMatrix::solveLy(vector<double>& b, vector<double>& y)
{
	y[0] = b[0];
	for (int i = 1; i < y.size(); i++)
	{
		y[i] = b[i] - y[i - 1] * bot[i];
	}
}

void BandMatrix::multAx(vector<double>& x)
{
	int size = x.size();
	double tmp[2];
	tmp[0] = diag[0] * x[0] + top[0] * x[1];
	tmp[1] = 0;
	for (int i = 1; i < size - 1; i++)
	{
		tmp[i%2] = bot[i] * x[i - 1] + diag[i] * x[i] + top[i] * x[i + 1];
		x[i - 1] = tmp[(i+1)%2];
	}
	tmp[size % 2] = bot[size] * x[size - 1] + diag[size] * x[size];
	x[size - 1] = tmp[(size + 1) % 2];
	x[size] = tmp[size % 2];
}