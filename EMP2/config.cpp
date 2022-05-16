#include "config.h"

const double eps = 1e-12;
const int iterMax = 1e5;
const double relaxParam = 0.7;

double getLambda(int i, double x, double u)
{
	return u;
}

double getGamma(int i, double x)
{
	return 1;
}

double getF(double x)
{
	return x - 1;
}

double getS1(double x)
{
	return x;
}