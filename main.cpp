#define BOOST_UBLAS_NDEBUG 1 
#include <iostream>
#include "source/marketdata.h"
#include "source/stylefactor.h"
#include "source/covariance.h"
using namespace std;
namespace ublas = boost::numeric::ublas;

int main(int argc, char* argv[])
{
	// retrive marketdata
	MarketData::instance().init_market_data();

	// construct a StyleFactor with given time range
	StyleFactor* SF = new StyleFactor(20070104, 20150104);

	// compute the barra style descriptors
	SF->compute_descriptors();
	// save the descriptors to csv
	SF->to_csv((*SF).size, "size_des");
	SF->to_csv((*SF).beta, "beta_des");
	SF->to_csv((*SF).hsigma, "hsigma_des");
	SF->to_csv((*SF).rstr, "rstr_des");
	SF->to_csv((*SF).dastd, "dastd_des");
	SF->to_csv((*SF).cmra, "cmra_des");
	SF->to_csv((*SF).stom, "stom_des");
	SF->to_csv((*SF).stoq, "stoq_des");
	SF->to_csv((*SF).stoa, "stoa_des");

	// compute factor exposures
	SF->compute_exposures();
	// save the factor exposures
	SF->to_csv((*SF).Size, "Size_exp");
	SF->to_csv((*SF).Beta, "Beta_exp");
	SF->to_csv((*SF).Momentum, "Momentum_exp");
	SF->to_csv((*SF).ResVol, "Res_Vol_exp");
	SF->to_csv((*SF).NlSize, "Nl_Size_exp");
	SF->to_csv((*SF).Liquidity, "Liquidity_exp");

	// compute factor returns
	SF->compute_returns();
	SF->rt_to_csv();
	
	// factor covariance matrix
	Covariance* COV = new Covariance(SF->factor_rt, SF->dates);
	//Covariance* COV = new Covariance(tmp_rt, tmp_d);
	ublas::matrix<double> vanilla_cov = COV->compute_vallina_cov();
	cout<<"Vanilla covariance matrix: "<<endl;
	cout<<vanilla_cov<<endl;
	ublas::matrix<double> nw_cov = COV->compute_nw_cov(5);
	cout<<"Newey-west covariance matrix: "<<endl;
	cout<<nw_cov<<endl;
	ublas::matrix<double> inverse(6, 6);
	COV->invert_cov(nw_cov, inverse);
	cout<<"Inversed covariance matrix: "<<endl;
	cout<<inverse<<endl;
	cout<<"Product of nw-cov and its inverse: "<<endl;
	cout<<ublas::prod(nw_cov, inverse)<<endl;

	return 1;
}

