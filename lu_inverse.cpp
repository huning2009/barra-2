// Created by Ruiqi Zhang on 2019-07-14

#include <iostream>
#include <vector>
using namespace std;


void print(const vector<vector<double>>& M)
{
	for (int i = 0; i != M.size(); ++i)
	{
		for (int j = 0; j != M[i].size(); ++j)
			cout<<M[i][j]<<"  ";
		cout<<endl;
	}
	cout<<"\n"<<endl;
}

// transpose of a matrix
vector<vector<double>> trans(const vector<vector<double>>& M)
{
	const int n = M.size();
	vector<vector<double>> T(n, vector<double>(n));
	for (int i = 0; i != n; ++i)
		for (int j = 0; j != n; ++j)
			T[i][j] = M[j][i];
	return T;
}

// product of square matrix
vector<vector<double>> prod(const vector<vector<double>>& A, vector<vector<double>>& B)
{
	const int n = A.size();
	vector<vector<double>> P(n, vector<double>(n));
	for (int i = 0; i != n; ++i)
	{
		for (int j = 0; j != n; ++j)
		{
			for (int k = 0; k != n; ++k)
				P[i][j] += A[i][k] * B[k][j];
		}
	}
	return P;
}

// inverse of a lower triangular matrix
vector<vector<double>> invertLT(const vector<vector<double>>& L)
{
	const int n = L.size();
	vector<vector<double>> M(n, vector<double>(n));
	for (int i = 0; i != n; ++i)
	{
		if (L[i][i] == 0.0)
		{
			cout<<"Input is singular matrix."<<endl;
			return M;
		}
		M[i][i] = 1.0 / L[i][i];
		for (int j = 0; j < i; ++j)
		{
			for (int k = j; k < i; ++k)
				M[i][j] += L[i][k] * M[k][j];
			M[i][j] = -M[i][j] / L[i][i];
		}
	}
	return M;
}

// LU decomposition 
void LU_Decompose(vector<vector<double>>& M, vector<vector<double>>& L, vector<vector<double>>& U)
{
	const int n = M.size();
	for (int i = 0; i != n; ++i)
	{
		// Upper triangular
		for (int k = i; k < n; ++k)
		{
			double sum = 0.;
			for (int j = 0; j < i; ++j)
				sum += (L[i][j] * U[j][k]);
			U[i][k] = M[i][k] - sum;
		}

		// Lower triangular
		for (int k = i; k < n; ++k)
		{
			if (i == k)
				L[i][i] = 1.;
			else
			{
				double sum = 0.;
				for (int j = 0; j < i; ++j)
					sum += (L[k][j] * U[j][i]);
				L[k][i] = (M[k][i] - sum) / U[i][i];
			}
		}
	}
	return;
}

vector<vector<double>> LU_Invert(vector<vector<double>>& M)
{
	const int n = M.size();
	vector<vector<double>> L(n, vector<double>(n));
	vector<vector<double>> U(n, vector<double>(n));
	LU_Decompose(M, L, U);
	vector<vector<double>> invL = invertLT(L);
	vector<vector<double>> invU = trans(invertLT(trans(U)));
	vector<vector<double>> invM = prod(invU, invL);
	return invM;
}

int main()
{
	// test cases for the functions
	// uncomment them if needed
	/*
	vector<vector<double>> L = {
		{1, 0, 0, 0},
		{4, 1, 0, 0},
		{2, 5, 1, 0},
		{6, 8, 4, 1}
	};
	
	vector<vector<double>> invL = invertLT(L);
	vector<vector<double>> PL = prod(L, invL);
	cout<<"Inverse of a lower triangular matrix: "<<endl;
	print(L);
	print(invL);
	print(PL);

	vector<vector<double>> U = {
		{1, 2, 1, 1},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}	
	};

	vector<vector<double>> invU = trans(invertLT(trans(U)));
	vector<vector<double>> PU = prod(U, invU);
	cout<<"Inverse of an upper triangular matrix: "<<endl;
	print(U);
	print(invU);
	print(PU);

	vector<vector<double>> M = {
		{2, -1, -2},
		{-4, 6, 3},
		{-4, -2, 8}
	};
	vector<vector<double>> ML(3, vector<double>(3));
	vector<vector<double>> MU(3, vector<double>(3));
	LU_Decompose(M, ML, MU);
	cout<<"LU decompose of test matrix: "<<endl;
	print(M);
	print(ML);
	print(MU);
	vector<vector<double>> PM = prod(ML, MU);
	print(PM);

	cout<<"Inverse matrix:"<<endl;
	print(M);
	vector<vector<double>> invM = LU_Invert(M);
	print(invM);
	PM = prod(M, invM);
	print(PM);
	*/

	vector<vector<double>> nw_cov = {
		{2.52434e-05,-5.129e-06,-8.94367e-06,-3.71576e-06,-4.57478e-06,-2.8802e-06},
		{-4.39538e-06,5.14014e-06,-2.44862e-06,1.38609e-06,3.90801e-07,-7.30288e-08},
		{-8.7635e-06,-2.8149e-06,1.68069e-05,-5.018e-06,-6.36429e-08,-1.46891e-07},
		{-3.0625e-06,1.74427e-06,-4.68974e-06,6.04223e-06,2.78349e-07,-3.12608e-07},
		{-5.39954e-06,7.18235e-07,-2.00015e-07,9.08732e-07,2.62369e-06,1.3489e-06},
		{-3.62248e-06,3.41258e-07,-5.24892e-07,3.96708e-07,1.34558e-06,2.06383e-06}
	};

	double lambda = 0.02;
	const int n = 6;
	vector<vector<double>> M(6, vector<double>(6));
	for (int i = 0; i != n; ++i)
	{
		for (int j = 0; j != n; ++j)
		{
			M[i][j] = 1000000. * nw_cov[i][j];
			if (i == j)
				M[i][j] += lambda;
		}
	}
	vector<vector<double>> invM = LU_Invert(M);
	for (int i = 0; i != n; ++i)
		for (int j = 0; j != n; ++j)
			invM[i][j] = 1000000. * invM[i][j];

	vector<vector<double>> CM = prod(M, invM);

	cout<<"Newey-west covariance matrix: "<<endl;
	print(nw_cov);
	cout<<"Inverse matrix of Newey-west covariance: "<<endl;
	print(invM);
	cout<<"Product: "<<endl;
	print(prod(nw_cov, invM));

}	
