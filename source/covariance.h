// Created by Ruiqi Zhang on 2019-07-11

#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include <map>

class Covariance
{
	public:
		std::vector<int> dates;
		std::vector<boost::numeric::ublas::matrix<double>> factors; // factor returns 
		boost::numeric::ublas::matrix<double> vanilla_cov; // vanilla factor covariance matrix
		boost::numeric::ublas::matrix<double> nw_cov; // factor covariance matrix with newey-west correction
		boost::numeric::ublas::matrix<double> inv_cov; // inverse matrix
			
	public:
		Covariance(std::vector<std::vector<double>>& vec_rt, std::vector<int>& ds);
		~Covariance();

		boost::numeric::ublas::matrix<double> compute_vallina_cov();
		boost::numeric::ublas::matrix<double> compute_nw_cov(int q);
		bool invert_cov(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);	
};


#endif // COVARIANCE_H
