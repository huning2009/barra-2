// Created by Ruiqi Zhang on 2019-07-08

#ifndef STYLEFACTOR_H
#define STYLEFACTOR_H

#include <map>
#include <vector>

class StyleFactor
{
	private:
		int begin_date;
		int end_date;
	public:
		StyleFactor(int b_date, int e_date);
		~StyleFactor();

	public:
		// descriptors
		std::map<int, std::map<std::string, double>> size;
		std::map<int, std::map<std::string, double>> beta;
		std::map<int, std::map<std::string, double>> hsigma;
		std::map<int, std::map<std::string, double>> rstr;
		std::map<int, std::map<std::string, double>> dastd;
		std::map<int, std::map<std::string, double>> cmra;
		std::map<int, std::map<std::string, double>> stom;
		std::map<int, std::map<std::string, double>> stoq;
		std::map<int, std::map<std::string, double>> stoa;

		// exposures
		std::map<int, std::map<std::string, double>> Size;
		std::map<int, std::map<std::string, double>> Beta;
		std::map<int, std::map<std::string, double>> Momentum;
		std::map<int, std::map<std::string, double>> ResVol;
		std::map<int, std::map<std::string, double>> NlSize;
		std::map<int, std::map<std::string, double>> Liquidity;

		// returns
		std::map<int, double> Size_Rt;
		std::map<int, double> Beta_Rt;
		std::map<int, double> Momentum_Rt;
		std::map<int, double> ResVol_Rt;
		std::map<int, double> NlSize_Rt;
		std::map<int, double> Liquidity_Rt;
		std::vector<std::vector<double>> factor_rt;
		std::vector<int> dates;

		void compute_descriptors();
		void compute_exposures();
		void compute_returns();

		bool fillna_mean(std::vector<double>& vec);
		bool fillna_mean(std::map<std::string, double>& input, std::vector<double>& cap_w);
		bool zscore(std::map<std::string, double>& input, std::vector<double>& cap_w);
		bool remove_extreme(std::map<std::string, double>& input);
		void normalize(std::map<std::string, double>& input, std::vector<double>& cap_w);
		void to_csv(std::map<int, std::map<std::string, double>>& input, std::string fname);
		void rt_to_csv(); // store factor returns
};

#endif // STYLEFACTOR_H
