#pragma once

extern const double eps;
extern const double delta;
extern const int iterMax;
extern const double relaxParam;

extern const bool isS2Left;
extern const bool isS2Right;
extern const bool isS3Left;
extern const bool isS3Right;

double getLambda(int i, double u,double x);

double getGamma(int i, double x);

double getF(int i, double x);

double getS1(int i, double x);

double getTheta(int i, double x, double u);

double getBeta(int i, double x, double u);

double getUb(int i, double x, double u);