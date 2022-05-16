#pragma once
#include <vector>

using namespace std;

class BandMatrix // 3 diagonals
{
public:
	vector<double> diag;
	vector<double> top;
	vector<double> bot;

	void resize(int size);
	void toZero();
	void matrixToLU(); //main diagonal of L consists of only 1
	void solveLUx(vector<double>& b, vector<double>& x);
	void multAx(vector<double>& x);

private:
	void solveUx(vector<double>& b, vector<double>& x);
	void solveLy(vector<double>& b, vector<double>& y);
};
