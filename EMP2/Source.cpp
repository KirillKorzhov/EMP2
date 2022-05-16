#include <iostream>
#include <fstream>
#include <locale>

#include "config.h"

using namespace std;

string const nodetxt = "node.txt";
string const fragmenttxt = "fragment.txt";
string const s1txt = "s1.txt";

class FEM
{
public:
	~FEM()
	{
		delete[] nodeArray;
		delete[] b;
		delete[] diag;
		delete[] top;
		delete[] bot;
		delete[] q;
		delete[] fragmentArray;
		delete[] s1Array;
	}

	void nodeRead()
	{
		ifstream in;
		in.open(nodetxt);
		in >> nodeCount;

		nodeArray = new double[nodeCount];
		b = new double[nodeCount];
		q = new double[nodeCount]();
		diag = new double[nodeCount];//у нас одномерна€ задача, так что матрица будет трехдиагональна€. ‘ормат матрица ленточный - это уже условие задачи
		top = new double[nodeCount];
		bot = new double[nodeCount]; 
		fragmentArray = new int[nodeCount];

		for (int i = 0; i < nodeCount; i++)
			in >> nodeArray[i];
		in.close();
	}

	void fragmentRead()
	{
		ifstream in;
		in.open(fragmenttxt);

		for (int i = 0; i < nodeCount; i++)
			in >> fragmentArray[i];

		in.close();
	}

	void s1Read()
	{
		ifstream in;
		in.open(s1txt);
		in >> s1Count;

		s1Array = new int[s1Count];

		for (int i = 0; i < s1Count; i++)
			in >> s1Array[i];

		in.close();
	}

	void iter()
	{
		globalToZero();
		makeGlobal();
		matrixToLU();
		solveLUx();
	}

	void printMatrix()
	{
		cout << endl;
		cout << fixed;
		cout.precision(2);

		cout << diag[0] << "  " << top[0];
		int tmp1 = 0;
		int tmp2 = nodeCount - 2;
		for (int i = 0; i < tmp2; i++)
		{
			cout << "  " << "0.00";
		}
		cout << endl;
		tmp2--;

		for (int i = 1; i < nodeCount - 1; i++)
		{
			for (int j = 0; j < tmp1; j++)
			{
				cout << "0.00" << "  ";
			}
			cout << bot[i] << "  " << diag[i] << "  " << top[i];
			for (int j = 0; j < tmp2; j++)
			{
				cout << "  " << "0.00";
			}
			cout << endl;
			tmp1++;
			tmp2--;
		}

		for (int j = 0; j < tmp1; j++)
		{
			cout << "0.00" << "  ";
		}
		cout << bot[nodeCount - 1] << "  " << diag[nodeCount - 1];
		cout << endl << endl;
		cout << defaultfloat;
	}

	void printQ()
	{
		for (int i = 0; i < nodeCount; i++)
		{
			cout << q[i] << "  ";
		}
		cout << endl;
	}


	void printB()
	{
		for (int i = 0; i < nodeCount; i++)
		{
			cout << b[i] << "  ";
		}
		cout << endl;
	}

private:
	int nodeCount = 0;
	double* nodeArray = NULL;
	double* b = NULL;
	double* diag = NULL;
	double* top = NULL;
	double* bot = NULL;
	double* q = NULL;
	int* fragmentArray = NULL;
	int s1Count = 0;
	int* s1Array = NULL;


	double getAvgLambda(int i)
	{
		return (getLambda(fragmentArray[i], nodeArray[i], q[i]) + getLambda(fragmentArray[i + 1], nodeArray[i + 1], q[i])) / 2;
	}

	double getAvgGamma(int i)
	{
		return (getGamma(fragmentArray[i], nodeArray[i]) + getGamma(fragmentArray[i + 1], nodeArray[i + 1])) / 2;
	}

	double getStep(int i)
	{
		return nodeArray[i + 1] - nodeArray[i];
	}
	
	void globalToZero()
	{
		vecToZero(b);
		vecToZero(diag);
		vecToZero(top);
		vecToZero(bot);
	}

	void vecToZero(void* vec)
	{
		memset(vec, 0, sizeof(vec) * nodeCount);
	}

	void makeGlobal()
	{
		makeLeft();
		makeRight();
		addS1();
	}

	void solveLUx()
	{
		solveLy(q, b);
		solveUx(q, q);

	}

	void matrixToLU()
	{
		for (int i = 1; i < nodeCount; i++)
		{
			bot[i] = bot[i] / diag[i - 1];
			diag[i] = diag[i] - bot[i] * top[i - 1];
		}
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
		for (int i = 0; i < nodeCount - 1; i++)
		{
			double tmp = getAvgLambda(i) / getStep(i);
			diag[i] += tmp;
			top[i] -= tmp;
			bot[i + 1] -= tmp;
			diag[i + 1] += tmp;
		}
	}

	void addMToLeft()
	{
		for (int i = 0; i < nodeCount - 1; i++)
		{
			double tmp = getAvgGamma(i) * getStep(i) / 6;
			diag[i] += 2 * tmp;
			top[i] += tmp;
			bot[i + 1] += tmp;
			diag[i + 1] += 2 * tmp;
		}
	}

	void addFToRight()
	{
		for (int i = 0; i < nodeCount - 1; i++)
		{
			b[i] += (2 * getF(nodeArray[i]) + getF(nodeArray[i + 1])) * getStep(i) / 6;
			b[i + 1] += (getF(nodeArray[i]) + 2 * getF(nodeArray[i + 1])) * getStep(i) / 6;
		}
	}

	void addS1()
	{
		for (int i = 0; i < s1Count; i++)
		{
			diag[s1Array[i]] = 1;
			top[s1Array[i]] = 0;
			bot[s1Array[i]] = 0;

			b[s1Array[i]] = getS1(nodeArray[s1Array[i]]);
		}
	}

	void solveUx(double* x, double* b)
	{
		x[nodeCount - 1] = b[nodeCount - 1];
		for (int i = nodeCount - 2; i >= 0; i--)
		{
			x[i] = (b[i] - top[i] * x[i + 1]) / diag[i];
		}
	}

	void solveLy(double* y, double* b)
	{
		y[0] = b[0];
		for (int i = 1; i < nodeCount; i++)
		{
			y[i] = b[i] - y[i - 1] * bot[i];
		}
	}

};

int main()
{
	setlocale(LC_ALL, "Russian");

	FEM fem;
	fem.nodeRead();
	fem.fragmentRead();
	fem.s1Read();

	for (int i = 0; i < 1000; i++)
	{
		fem.iter();
	}

	fem.printQ();

	return 0;
}