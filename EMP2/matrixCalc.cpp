#include "matrixCalc.h"
#include <vector>

using namespace std;

double norm(vector<double>& x)
{
	double tmp = 0;
	for (int i = 0; i < x.size(); i++)
	{
		tmp += x[i] * x[i];
	}
	return sqrt(tmp);
}

void alphaMULTx(const double alpha, vector<double>& x)
{
	for (int i = 0; i < x.size(); i++)
	{
		x[i] *= alpha;
	}
}

void vecADDvec(vector<double>& x, vector<double>& y)
{
	for (int i = 0; i < x.size(); i++)
	{
		x[i] += y[i];
	}
}