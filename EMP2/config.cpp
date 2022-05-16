#include <iostream>
#include "config.h"

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
	return x*x - 6*x*x;
}

double getS1(double x)
{
	return x*x;
}