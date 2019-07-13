// Created by Ruiqi Zhang on 2019-07-11
#define BOOST_UBLAS_NDEBUG 1 
#include "covariance.h"
#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace std;
namespace ublas = boost::numeric::ublas;

Covariance::Covariance(vector<vector<double>>& vec_rt, vector<int>& ds)
{
	dates.clear();
	factors.clear();
	// init member variables
	for (int i = 0; i != vec_rt.size(); ++i)
	{
		int date = ds[i];
		ublas::matrix<double> tmp_matrix(6, 1);
		tmp_matrix.clear();
		// demean 
		double sum = 0.;
		for (int j = 0; j != vec_rt[i].size(); ++j)
			sum += vec_rt[i][j];
		for (int j = 0; j != vec_rt[i].size(); ++j)
			tmp_matrix(j, 0) = (vec_rt[i][j] - sum / double(vec_rt[i].size()));
		dates.push_back(date);
		factors.push_back(tmp_matrix);
	} // loop of dates
}

Covariance::~Covariance() {};

ublas::matrix<double> Covariance::compute_vallina_cov() 
{
	ublas::matrix<double> Vf(6, 6);
	Vf.clear();
	for (int t = 0; t != factors.size(); ++t)
	{	
		ublas::matrix<double> Ft = factors[t];
		Vf += ublas::prod(Ft, trans(Ft));
	}
	cout<<"Vf: "<<endl;
	cout<<Vf<<endl;
	Vf /= double(factors.size());
	return Vf;
};

ublas::matrix<double> Covariance::compute_nw_cov(int q)
{
	ublas::matrix<double> Vf(6, 6);
	Vf.clear();
	for (int i = 0; i <= q; ++i)
	{
		for (int  t = (factors.size() - 1); t >= i; --t)
		{
			ublas::matrix<double> Ft = factors[t];
			ublas::matrix<double> Fti = factors[t - i];
			double w = 1. - (i / double(1 + q));
			Vf += w * ublas::prod(Ft, trans(Fti)) / double(factors.size());
		}
	} // loop of lags
	return Vf;
}

bool Covariance::invert_cov(ublas::matrix<double>& input, ublas::matrix<double>& inverse)
{
	ublas::matrix<double> A(input);
	// avoid singularity
	A = 1000000. * A;
	double lambda = 1.;
	for (int i = 0; i != A.size1(); ++i)
		A(i, i) += lambda;
	ublas::permutation_matrix<std::size_t> pm(A.size1());
	// perform L-U factorization
	int res = lu_factorize(A,pm);
	if (res != 0)
		return false;
	// create identity matrix of "inverse"
    inverse.assign(ublas::identity_matrix<double>(A.size1()));
	// backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);
	// 
	inverse = inverse * 1000000.;
	return true;
}




