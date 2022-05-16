#include <iostream>
#include <fstream>
#include <locale>
#include <vector>

#include "config.h"
#include "BandMatrix.h"
#include "matrixCalc.h"

using namespace std;

string const nodetxt = "node.txt";
string const fragmenttxt = "fragment.txt";
string const s1txt = "s1.txt";

class FEM
{
public:
	void nodeRead()
	{
		ifstream in;
		in.open(nodetxt);
		int nodeCount;
		in >> nodeCount;

		nodeVec.resize(nodeCount);
		b.resize(nodeCount);
		q.resize(nodeCount);
		qPrev.resize(nodeCount);
		bandMatrix.resize(nodeCount);
		fragmentVec.resize(nodeCount);

		for (int i = 0; i < nodeCount; i++)
			in >> nodeVec[i];
		in.close();
	}

	void fragmentRead()
	{
		ifstream in;
		in.open(fragmenttxt);

		for (int i = 0; i < nodeVec.size(); i++)
			in >> fragmentVec[i];

		in.close();
	}

	void s1Read()
	{
		ifstream in;
		in.open(s1txt);
		int s1Count;
		in >> s1Count;

		s1Vec.resize(s1Count);

		for (int i = 0; i < s1Count; i++)
			in >> s1Vec[i];

		in.close();
	}

	int simpleIter()
	{
		addFToRight();
		double normDenominator = norm(b);
		int i = 0;
		for (; i < iterMax; i++)
		{
			double dis = iter(normDenominator);
			alphaMULTx(relaxParam, q);
			alphaMULTx(1 - relaxParam, qPrev);
			vecADDvec(q, qPrev);
			if (dis < eps)
			{
				break;
			}
			cout << dis << endl;
			q.swap(qPrev);
		}
		return i;
	}

	double iter(double normDenominator)
	{
		globalToZero();
		makeGlobal();
		double prevDiscrepancy = discrepancy(normDenominator);
		bandMatrix.matrixToLU();
		bandMatrix.solveLUx(b, q);
		return prevDiscrepancy;
	}

	void printQ()
	{
		for (int i = 0; i < q.size(); i++)
		{
			cout << q[i] << "  ";
		}
		cout << endl;
	}

private:
	vector <double> nodeVec;
	BandMatrix bandMatrix;
	vector <double>	b;
	vector <double>	q;
	vector <double> qPrev;
	vector <int> fragmentVec;
	vector <int> s1Vec;

	double getAvgLambda(int i)
	{
		return (getLambda(fragmentVec[i], nodeVec[i], qPrev[i]) + getLambda(fragmentVec[i + 1], nodeVec[i + 1], qPrev[i + 1])) / 2;
	}

	double getAvgGamma(int i)
	{
		return (getGamma(fragmentVec[i], nodeVec[i]) + getGamma(fragmentVec[i + 1], nodeVec[i + 1])) / 2;
	}

	double getStep(int i)
	{
		return nodeVec[i + 1] - nodeVec[i];
	}
	
	void globalToZero()
	{
		fill(b.begin(), b.end(), 0);
		bandMatrix.toZero();
	}

	void makeGlobal()
	{
		makeLeft();
		makeRight();
		addS1();
	}

	void makeLeft()
	{
		addGToLeft();
		addMToLeft();
	}

	void makeRight()
	{
		addFToRight();
	}

	void addGToLeft()
	{
		for (int i = 0; i < nodeVec.size() - 1; i++)
		{
			double tmp = getAvgLambda(i) / getStep(i);
			bandMatrix.diag[i] += tmp;
			bandMatrix.top[i] -= tmp;
			bandMatrix.bot[i + 1] -= tmp;
			bandMatrix.diag[i + 1] += tmp;
		}
	}

	void addMToLeft()
	{
		for (int i = 0; i < nodeVec.size() - 1; i++)
		{
			double tmp = getAvgGamma(i) * getStep(i) / 6;
			bandMatrix.diag[i] += 2 * tmp;
			bandMatrix.top[i] += tmp;
			bandMatrix.bot[i + 1] += tmp;
			bandMatrix.diag[i + 1] += 2 * tmp;
		}
	}

	void addFToRight()
	{
		for (int i = 0; i < nodeVec.size() - 1; i++)
		{
			b[i] += (2 * getF(nodeVec[i]) + getF(nodeVec[i + 1])) * getStep(i) / 6;
			b[i + 1] += (getF(nodeVec[i]) + 2 * getF(nodeVec[i + 1])) * getStep(i) / 6;
		}
	}

	void addS1()
	{
		for (int i = 0; i < s1Vec.size(); i++)
		{
			bandMatrix.diag[s1Vec[i]] = 1;
			bandMatrix.top[s1Vec[i]] = 0;
			bandMatrix.bot[s1Vec[i]] = 0;

			b[s1Vec[i]] = getS1(nodeVec[s1Vec[i]]);
		}
	}

	double discrepancy()
	{
		return normAqSUBb() / norm(b);
	}

	double discrepancy(double normDenominator)
	{
		return normAqSUBb() / normDenominator;
	}

	double normAqSUBb()
	{
		int end = qPrev.size() - 1;
		double tmp = (bandMatrix.diag[0] * qPrev[0] + bandMatrix.top[0] * qPrev[1] - b[0]);
		double result = tmp * tmp;
		for (int i = 1; i < end; i++)
		{
			tmp = bandMatrix.bot[i] * qPrev[i - 1] + bandMatrix.diag[i] * qPrev[i] + bandMatrix.top[i] * qPrev[i + 1] - b[i];
			result += tmp * tmp;
		}
		tmp = bandMatrix.bot[end] * qPrev[end - 1] + bandMatrix.diag[end] * qPrev[end] - b[end];
		result += tmp * tmp;
		return sqrt(result);
	}

};

int main()
{
	setlocale(LC_ALL, "Russian");

	FEM fem;
	fem.nodeRead();
	fem.fragmentRead();
	fem.s1Read();

	cout << fem.simpleIter() << endl;

	fem.printQ();

	return 0;
}