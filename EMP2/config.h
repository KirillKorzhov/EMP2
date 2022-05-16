#pragma once

extern const double eps;
extern const int iterMax;
extern const double relaxParam;

double getLambda(int i, double u,double x);

double getGamma(int i, double x);

double getF(double x);

double getS1(double x);