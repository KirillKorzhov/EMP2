#include "config.h"

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
	return x*x-2;
}

double getS1(int i, double x)
{
	return x*x;
}

double getTheta(int i, double x, double u)
{
	return 2*x*getLambda(0,x,u);
}

double getBeta(int i, double x, double u)
{
	return 1;
}

double getUb(int i, double x, double u)
{
	return 1;
}
