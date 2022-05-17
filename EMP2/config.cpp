#include "config.h"
#include <math.h>
const double eps = 1e-12;
const double delta = 1e-12;
const int iterMax = 1e4;
const double relaxParam = 1;

const bool isS2Left = false;
const bool isS2Right = true;
const bool isS3Left = false;
const bool isS3Right = false;


double getLambda(int i, double x, double u)
{
	return 1;
}

double getGamma(int i, double x)
{
	return 1;
}

double getF(int i, double x)
{
	return x;
}

double getS1(int i, double x)
{
	return x;
}

double getTheta(int i, double x, double u)
{
	return 1;
}

double getBeta(int i, double x, double u)
{
	return 1;
}

double getUb(int i, double x, double u)
{
	return 1;
}
