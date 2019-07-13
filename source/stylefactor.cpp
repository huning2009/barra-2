// Created by Ruiqi Zhang on 2019-07-08

#include "stylefactor.h"
#include "marketdata.h"
#include "gsl/gsl_fit.h"
#include "gsl/gsl_multifit.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
using namespace std;

StyleFactor::StyleFactor(int b_date, int e_date)
{
	begin_date = b_date;
	end_date = e_date;

	size.clear();
	beta.clear();
	hsigma.clear();
	rstr.clear();
	dastd.clear();
	cmra.clear();
	stom.clear();
	stoq.clear();
	stoa.clear();

	Size.clear();
	Beta.clear();
	Momentum.clear();
	ResVol.clear();
	NlSize.clear();
	Liquidity.clear();

	Size_Rt.clear();
	Beta_Rt.clear();
	Momentum_Rt.clear();
	ResVol_Rt.clear();
	NlSize_Rt.clear();
	Liquidity_Rt.clear();
	factor_rt.clear();
	dates.clear();
}

StyleFactor::~StyleFactor() {};

void StyleFactor::compute_descriptors()
{
	auto start = chrono::high_resolution_clock::now();
	// retrieve data from MarketData::instance()
	int N = MarketData::instance().trading_days.size();
	int S = MarketData::instance().symbols.size();
	for (int i = 0; i != N; ++i)
	{
		if (MarketData::instance().trading_days[i] < begin_date || MarketData::instance().trading_days[i] > end_date)
			continue;
		int date = MarketData::instance().trading_days[i];
		// set the structure of the maps
		map<string, double> temp_map;
		temp_map.clear();
		for (int s = 0; s != S; ++s)
			temp_map[MarketData::instance().symbols[s]] = numeric_limits<double>::quiet_NaN();
		// size
		map<string, double> size_on_date = temp_map;
		//size_on_date.clear();
		// beta
		map<string, double> beta_on_date = temp_map;
		//beta_on_date.clear(); 
		// hsigma
		map<string, double> hsigma_on_date = temp_map;
		//hsigma_on_date.clear();
		// rstr
		map<string, double> rstr_on_date = temp_map;
		//rstr_on_date.clear();
		// dastd
		map<string, double> dastd_on_date = temp_map;
		//dastd_on_date.clear();
		// cmra 
		map<string, double> cmra_on_date = temp_map;
		//cmra_on_date.clear();
		// stom
		map<string, double> stom_on_date = temp_map;
		//stom_on_date.clear();
		// stoq
		map<string, double> stoq_on_date = temp_map;
		//stoq_on_date.clear();
		// stoa
		map<string, double> stoa_on_date = temp_map;
		//stoa_on_date.clear();	
		// loop of symbols
		for (int s = 0; s != S; ++s)
		{
			string symbol = MarketData::instance().symbols[s];
			//size
			if (isnormal(MarketData::instance().market_data[date].tot_mkt_val[s]))
				size_on_date[symbol] = log(MarketData::instance().market_data[date].tot_mkt_val[s]);
			if (i < (N - 252))
			{
				// beta
				vector<double> X, Y, W;
				X.clear(); Y.clear(); W.clear();
				// dastd
				vector<double> daily_excess_rt;
				daily_excess_rt.clear();
				// cmra
				vector<double> vec_rt, vec_rft;
				vec_rt.clear(); vec_rft.clear();
				// stom, stoq, stoa
				vector<double> vec_share, vec_vol, vec_turnover;
				vec_share.clear(); vec_vol.clear(); vec_turnover.clear();
				// flags for nan values
				bool rt_na = false;
				bool rft_na = false;
				bool shr_na = false;
				bool vol_na = false;
				bool to_na = false;
				// loop for the time window
				for (int d = 251; d >= 0; --d)	
				{
					int f_date = MarketData::instance().trading_days[i + d];

					double rt = MarketData::instance().market_data[f_date].returns[s];
					double rft = MarketData::instance().market_data[f_date].shibor;
					double shr = MarketData::instance().market_data[f_date].shares[s];
					double vol = MarketData::instance().market_data[f_date].volume[s];
					double turnover = MarketData::instance().market_data[f_date].turnover[s];
					if (!isnormal(rt))
						rt_na = true;
					if (!isnormal(rft))
						rft_na = true;
					if (!isnormal(shr))
						shr_na = true;
					if (!isnormal(vol))
						vol_na = true;
					if (!isnormal(turnover))
						to_na = true;

					if (rft_na)
						cout<<"Found nan in Shibor!"<<endl;
					X.push_back(MarketData::instance().market_data[f_date].mkt_return - rft);
					W.push_back(pow(pow(2, -1/63), 252 - d));	

					vec_rt.push_back(rt);
					vec_rft.push_back(rft);
					
					vec_share.push_back(shr);
					vec_vol.push_back(vol);
					vec_turnover.push_back(turnover);
				}

				// process the nan values
				if (rt_na)
					rt_na = !fillna_mean(vec_rt);
				if (rft_na)
					rft_na = !fillna_mean(vec_rft);
				if (shr_na)
					shr_na = !fillna_mean(vec_share);
				if (vol_na)
					vol_na = !fillna_mean(vec_vol);
				if (to_na)
					to_na = !fillna_mean(vec_turnover);

				for (int d = 0; d != 252; ++d)
				{
					Y.push_back(vec_rt[d] - vec_rft[d]);
					double w_dastd = pow(pow(2, -1/42), d + 1);
                    daily_excess_rt.push_back(w_dastd * (vec_rt[d] - vec_rft[d]));
				}

				// WLS fit for beta and hsigma 
				if (X.size() == 252 && Y.size() == 252 && W.size() == 252 && !rt_na && !rft_na)
				{
					double c0, c1, cov00, cov01, cov11, chisq;
					gsl_fit_wlinear(&X[0], 1, &W[0], 1, &Y[0], 1, 252, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
					beta_on_date[symbol] = c1;
					// residuals for hsigma
					double resrt_sum = 0;
                    double resrt_square_sum = 0.;
					for (int t = 0; t != 252; ++t)
					{
						double resrt = Y[t] - (c0 + X[t] * c1);
                        resrt_sum += W[t] * resrt;
                        resrt_square_sum += W[t] * resrt * W[t] * resrt;
					}
					double resrt_var = resrt_square_sum / 252. - resrt_sum * resrt_sum / (252. * 252.);
					if (resrt_var < 0)
                        cout<<"Variance of resrt is negative."<<endl;
					else 
						hsigma_on_date[symbol] = sqrt(resrt_var);
				}

				// dastd
				if (daily_excess_rt.size() == 252 && !rt_na)
				{
					double sum = 0;
	                double sum_square = 0;
    	            for (auto r : daily_excess_rt)
        	        {
            	        sum += r;
                	    sum_square += r * r;
                	}
                	double var = sum_square / 252. - sum * sum / (252. * 252.);
					if (var < 0)
					{
	                    cout<<"Variance of daily excess is negative."<<endl;
						dastd_on_date[symbol] = numeric_limits<double>::quiet_NaN();
					}
					else
						dastd_on_date[symbol] = sqrt(var);
				}

				// cmra
				if (vec_rt.size() == 252 && vec_rft.size() == 252 && !rt_na && !rft_na)
				{
					// compounded monthly returns
					vector<double> cmp_rt, cmp_rft;
	                cmp_rt.clear();
    	            cmp_rft.clear();
					double cum_rt = 1., cum_rft = 1.;
					for (int d = 0; d != 252; d++)
                	{
                    	cum_rt *= (1 + vec_rt[d]);
	                    cum_rft *= (1 + vec_rft[d]);
    	                if ((d + 1) % 21 == 0)
        	            {
            	            cmp_rt.push_back(cum_rt - 1);
                	        cmp_rft.push_back(cum_rft - 1);
                        	cum_rt = 1; cum_rft = 1;
                    	}
                	}
					if (cmp_rt.size() != 12 || cmp_rft.size() != 12)
						cout<<"Size of cmp_rt or cmp_rft is not 12."<<endl;
					double Zmax = -9999.;
	                double Zmin = 9999.;
    	            double zt = 0;
        	        for (int t = 0; t != 12; ++t)
            	    {
                	    zt += (log(1 + cmp_rt[t]) - log(1 + cmp_rft[t]));
                    	if (zt > Zmax) Zmax = zt;
	                    if (zt < Zmin) Zmin = zt;
    	            }
					cmra_on_date[symbol] = log(1 + Zmax) - log(1 + Zmin);
				}

				// stom, stoq, stoa
				if (vec_turnover.size() == 252 && !to_na)
				{
					vector<double> vec_stom;
                    vec_stom.clear();
                    double tmp_sum = 0;
                    for (int t = 0; t != 252; ++t)
                    {
						tmp_sum += vec_turnover[t];
						if ((t + 1) % 21 == 0)
						{
							vec_stom.push_back(tmp_sum);
							tmp_sum = 0;
						}
					}
					stom_on_date[symbol] = log(vec_stom[0]);
                    int T = 3;
                    tmp_sum = 0;
                    for (int t = 0; t != T; ++t)
                        tmp_sum += (vec_stom[t]);
                    stoq_on_date[symbol] = log(tmp_sum / double(T));
                    T = 12;
                    tmp_sum = 0;
                    for (int t = 0; t != T; ++t)
                        tmp_sum += (vec_stom[t]);
                    stoa_on_date[symbol] = log(tmp_sum / double(T));
				}
				/*
				if (vec_share.size() == 252 && vec_vol.size() == 252 && !shr_na && !vol_na)
				{
					vector<double> vec_stom;
					vec_stom.clear();
					double tmp_sum = 0;
					for (int t = 0; t != 252; ++t)
					{
						if ((t + 1) % 21 == 0) 
						{
							tmp_sum += (vec_vol[t] / vec_share[t]);
							vec_stom.push_back(tmp_sum);
							tmp_sum = 0;
							continue;
						}
						tmp_sum += (vec_vol[t] / vec_share[t]);
					}
					stom_on_date[symbol] = log(vec_stom[0]);
					int T = 3;
					tmp_sum = 0;
					for (int t = 0; t != T; ++t)
						tmp_sum += (vec_stom[t]);
					stoq_on_date[symbol] = log(tmp_sum / double(T));
					T = 12;
					tmp_sum = 0;
					for (int t = 0; t != T; ++t)
                        tmp_sum += (vec_stom[t]);
                    stoa_on_date[symbol] = log(tmp_sum / double(T));
				}
				*/
			} // trailing 252 days

			// rstr
			if (i < (N - 504 - 21))
			{
				bool rt_na = false;
				bool rft_na = false;
				vector<double> vec_rt, vec_rft, vec_weight;
				vec_rt.clear(); vec_rft.clear(); vec_weight.clear();
				for (int d = 503; d >= 0; --d)
				{
					int f_date = MarketData::instance().trading_days[i + d];
					double rt = MarketData::instance().market_data[f_date].returns[s];
                    double rft = MarketData::instance().market_data[f_date].shibor;
					if (!isnormal(rt))
						rt_na = true;
					if (!isnormal(rft))
						rft_na = true;
					//double weight = pow(pow(2, -1/126), d - 21 + 1);
					//vec_sum.push_back(weight * (log(1 + rt) - log(1 + rft)));
					vec_weight.push_back(pow(pow(2, -1/126), 504 - d));
					vec_rt.push_back(rt);
					vec_rft.push_back(rft);
				}
				if (rt_na)
					rt_na = !fillna_mean(vec_rt);
				if (rft_na)
					rft_na = !fillna_mean(vec_rt);

				//vector<double> vec_sum;
                //vec_sum.clear();
				//for (int d = 0; d != 504; ++d)
				//	vec_sum.push_back(vec_weight[d] * (log(1 + vec_rt[d]) - log(1 + vec_rft[d])));

				if (vec_rt.size() == 504 && vec_rft.size() == 504 && vec_weight.size() == 504 && !rt_na && !rft_na)
				{
					double sum = 0.;
					for (int d = 0; d != 504; ++d)
						sum += vec_weight[d] * (log(1 + vec_rt[d]) - log(1 + vec_rft[d]));
					rstr_on_date[symbol] = sum;
				}
			} // trailing 504 + 21 days
		} // loop of symbols
		normalize(size_on_date, MarketData::instance().market_data[date].cap_weight);
		size[MarketData::instance().trading_days[i]] = size_on_date;
		if (i < (N - 252))
		{
			// normalize the descriptors
			normalize(beta_on_date, MarketData::instance().market_data[date].cap_weight);
			normalize(hsigma_on_date, MarketData::instance().market_data[date].cap_weight);
			normalize(dastd_on_date, MarketData::instance().market_data[date].cap_weight);
			normalize(cmra_on_date, MarketData::instance().market_data[date].cap_weight);
			normalize(stom_on_date, MarketData::instance().market_data[date].cap_weight);
			normalize(stoq_on_date, MarketData::instance().market_data[date].cap_weight);
			normalize(stoa_on_date, MarketData::instance().market_data[date].cap_weight);

			beta[MarketData::instance().trading_days[i + 251]] = beta_on_date;
			hsigma[MarketData::instance().trading_days[i + 251]] = hsigma_on_date;
			dastd[MarketData::instance().trading_days[i + 251]] = dastd_on_date;
			cmra[MarketData::instance().trading_days[i + 251]] = cmra_on_date;
			stom[MarketData::instance().trading_days[i + 251]] = stom_on_date;
			stoq[MarketData::instance().trading_days[i + 251]] = stoq_on_date;
			stoa[MarketData::instance().trading_days[i + 251]] = stoa_on_date;
		}
		if (i < (N - 504 - 21))
		{
			normalize(rstr_on_date, MarketData::instance().market_data[date].cap_weight);
			rstr[MarketData::instance().trading_days[i + 504 + 20]] = rstr_on_date;
		}

		//if ((i + 1) % 100 == 0)
		//	cout<<(i + 1)<<" days processed."<<endl;
	} // loop of trading days
	auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Time taken by compute_descriptors: "<<(duration.count() / 1000000.)<<" seconds"<<endl;
}

// vanilla mean
bool StyleFactor::fillna_mean(vector<double>& vec)
{
	if (vec.size() == 0)
		return false;
	int n = 0;
	double sum = 0.;
	for (int i = 0; i != vec.size(); ++i)
	{
		if (isnormal(vec[i]))
		{
			n++;
			sum += vec[i];
		}
	}
	// at least 70% valid values required 
	if (n / double(vec.size()) < 0.7) 
		return false;
	double mean = sum / double(n);
	for (int i = 0; i != vec.size(); ++i)
		if (!isnormal(vec[i]))
			vec[i] = mean;
	return true;
}

// cap-weighted mean
bool StyleFactor::fillna_mean(map<std::string, double>& input, vector<double>& cap_w)
{
	if (input.size() != cap_w.size())
	{
		cout<<"Inequal size between Input and cap_w."<<endl;
		return false;
	}

	// cap-weight mean
	double mean = 0.;
	int n = 0;
	int j = 0;
	for (auto i = input.begin(); i != input.end(); ++i)
	{
		if (isnormal(i->second))
		{
			mean += cap_w[j] * i->second;
 			n++;	   
		}
		j++;
	}
	
	//if (input.size() == 0 || n / double(input.size()) < 0.7)
	//	return false;
	
	for (auto i = input.begin(); i != input.end(); ++i)
//		i->second = 111.;
		if (!isnormal(i->second))
			i->second = mean;
	return true;
}

void StyleFactor::to_csv(map<int, map<string, double>>& input, string output)
{
	ofstream f_out("/Users/ruiqizhang/workarea/barra/output/" + output + ".csv");
	// feed column names
	f_out<<",";
	for (int i = 0; i != MarketData::instance().symbols.size(); ++i)
	{
		if ((i + 1) == MarketData::instance().symbols.size())
			f_out<<MarketData::instance().symbols[i];
		else
			f_out<<MarketData::instance().symbols[i]<<",";
	}
	f_out<<endl;

	for (int i = 0; i != MarketData::instance().trading_days.size(); ++i)
	{
		int date = MarketData::instance().trading_days[i];
		f_out<<MarketData::instance().trading_days[i]<<",";
		for (int j = 0; j != MarketData::instance().symbols.size(); ++j)
		{
			string symbol = MarketData::instance().symbols[j];
			if ((j + 1) == MarketData::instance().symbols.size())
			{
				if (input[date].find(symbol) == input[date].end())
					f_out<<""<<endl;
				else
					f_out<<input[date][symbol]<<endl;
			}
			else
			{
				if (input[date].find(symbol) == input[date].end())
					f_out<<",";
				else
					f_out<<input[date][symbol]<<",";
			}
		}
	}
	f_out.close();
}

void StyleFactor::rt_to_csv()
{
	ofstream f_out("/Users/ruiqizhang/workarea/barra/output/factor_rt.csv");
	// set column names
	f_out<<",Size,Beta,Momentum,ResVol,NlSize,Liquidity"<<endl;
	for (auto i = Size_Rt.begin(); i != Size_Rt.end(); ++i)
	{
		int date = i->first;
		f_out<<date<<","<<i->second<<","<<Beta_Rt[date]<<","<<Momentum_Rt[date]<<","<<ResVol_Rt[date]<<","<<NlSize_Rt[date]<<","<<Liquidity_Rt[date]<<endl;
	}
	f_out.close();
}

bool StyleFactor::zscore(map<string, double>& input, vector<double>& cap_w)
{
	if (input.size() != cap_w.size())
    {
        cout<<"Inequal size between input and cap weight: "<<input.size()<<"  "<<cap_w.size()<<endl;
        return false;
    }

    // cap-weight mean and normal std
    double mean = 0.;
    double sum = 0., square_sum = 0.;
    int n = input.size();
	int j = 0;
    for (auto i = input.begin(); i != input.end(); ++i)
    {
        if (isnormal(i->second))
        {
            mean += cap_w[j] * i->second;
            sum += i->second;
            square_sum += i->second * i->second;
        }
		++j;
    }
    if (sum < 0.00001)
        return false;
    
	double std = sqrt(square_sum / double(n) - sum * sum / double(n * n));
	if (!isnormal(std) || std < 0.000001)
	{
		cout<<"Abnormal std value: "<<std<<endl;
		return false;
	}

	// zscore
	for (auto i = input.begin(); i != input.end(); ++i)
	{
		if (isnormal(i->second))
			i->second = (i->second - mean) / std;
	}
	return true;
}

bool StyleFactor::remove_extreme(map<string, double>& input)
{
	vector<double> X;
    X.clear();
    for (auto j = input.begin(); j != input.end(); ++j)
    {
        if (isnormal(j->second))
            X.push_back(j->second);
    }
    if (X.size() == 0)
        return false;
    sort(X.begin(), X.end());
    double Xm = (X.size() % 2 == 0) ? 0.5 * (X[X.size() / 2 - 1] + X[X.size() / 2]) : X[X.size() / 2];
	if (!isnormal(Xm))
	{
		cout<<"Abnormal Xm value: "<<Xm<<endl;
		return false;
	}
    vector<double> Xmad;
    Xmad.clear();
    for (int s = 0; s != X.size(); ++s)
        Xmad.push_back(fabs(X[s] - Xm));
	sort(Xmad.begin(), Xmad.end());
    double Dmad = (Xmad.size() % 2 == 0) ? 0.5 * (Xmad[Xmad.size() / 2 - 1] + Xmad[Xmad.size() / 2]) : Xmad[Xmad.size() / 2];
	if (!isnormal(Dmad))
	{
		cout<<"Abnormal value Dmad value: "<<Dmad<<" "<<Xmad.size()<<" "<<Xmad[Xmad.size() / 2 - 1]<<" "<<Xmad[Xmad.size() / 2]<<endl;
		return false;
	}
    // deal with the extreme values
    double sig = 3.5;
	for (auto j = input.begin(); j != input.end(); ++j)
    {
        if (isnormal(j->second))
        {
            if (j->second > (Xm + sig * Dmad))
                j->second = (Xm + sig * Dmad);
            else if (j->second < (Xm - sig * Dmad))
                j->second = (Xm - sig * Dmad);
        }
    }
	return true;
}

void StyleFactor::normalize(map<string, double>& input, vector<double>& cap_w)
{
	zscore(input, cap_w);
	remove_extreme(input);
	zscore(input, cap_w);
}

void StyleFactor::compute_exposures()
{
	auto start = chrono::high_resolution_clock::now();
	// start from 504 + 21
	for (int i = (504 + 21); i != MarketData::instance().trading_days.size(); ++i)
	{
		int date = MarketData::instance().trading_days[i];
		// set the structure of the maps
        map<string, double> temp_map;
        temp_map.clear();
        for (int s = 0; s != MarketData::instance().symbols.size(); ++s)
            temp_map[MarketData::instance().symbols[s]] = numeric_limits<double>::quiet_NaN();
		// Size
		Size[date] = temp_map;
		// Beta
		Beta[date] = temp_map;
		// Momentum
		Momentum[date] = temp_map;
		// ResVol
		ResVol[date] = temp_map;
		// NlSize
		NlSize[date] = temp_map;
		// Liquidity
		Liquidity[date] = temp_map;

		// vectors to store size values
		vector<double> vec_size, vec_cubsize;
		vec_size.clear(); vec_cubsize.clear();
		vector<string> vec_tmp;
		vec_tmp.clear();
		for (int s = 0; s != MarketData::instance().symbols.size(); ++s)
		{
			string symbol = MarketData::instance().symbols[s];
			// retrieve size
			if (size.find(date) != size.end() && size[date].find(symbol) != size[date].end()) 
			{
				if (isnormal(size[date][symbol]))
				{
					vec_size.push_back(size[date][symbol]);
					vec_cubsize.push_back(pow(size[date][symbol], 3));
					vec_tmp.push_back(symbol);
					// fill Size factor exposure
					Size[date][symbol] = size[date][symbol];
				}
			}
			// retrieve beta
			if (beta.find(date) != beta.end() && beta[date].find(symbol) != beta[date].end())
			{
				if (isnormal(beta[date][symbol]))
					// fill Beta factor exposure
					Beta[date][symbol] = beta[date][symbol];
			}
			// retrieve rstr
			if (rstr.find(date) != rstr.end() && rstr[date].find(symbol) != rstr[date].end())
			{
				if (isnormal(rstr[date][symbol]))
					// fill Momentum factor exposure
					Momentum[date][symbol] = rstr[date][symbol];
			}
			// retrieve dastd, cmra and hsigma
			if (dastd.find(date) != dastd.end() && dastd[date].find(symbol) != dastd[date].end() && cmra.find(date) != cmra.end() && cmra[date].find(symbol) != cmra[date].end() && hsigma.find(date) != hsigma.end() && hsigma[date].find(symbol) != hsigma[date].end())
			{
				if (!isnormal(dastd[date][symbol]) && !isnormal(cmra[date][symbol]) && !isnormal(hsigma[date][symbol]))
					ResVol[date][symbol] = numeric_limits<double>::quiet_NaN();
				else
				{
					// fill ResVol factor exposure
					double tmp_dastd = isnormal(dastd[date][symbol]) ? dastd[date][symbol] : 0.;
					double tmp_cmra = isnormal(cmra[date][symbol]) ? cmra[date][symbol] : 0.;
					double tmp_hsigma = isnormal(hsigma[date][symbol]) ? hsigma[date][symbol] : 0.;
					ResVol[date][symbol] = 0.74 * tmp_dastd + 0.16 * tmp_cmra + 0.1 * tmp_hsigma;
				}
			}
			// retrieve stom, stoq, stoa
			if (stom.find(date) != stom.end() && stom[date].find(symbol) != stom[date].end() && stoq.find(date) != stoq.end() && stoq[date].find(symbol) != stoq[date].end() && stoa.find(date) != stoa.end() && stoa[date].find(symbol) != stoa[date].end())
			{
				if (!isnormal(stom[date][symbol]) && !isnormal(stoq[date][symbol]) && !isnormal(stoa[date][symbol]))
					Liquidity[date][symbol] = numeric_limits<double>::quiet_NaN();
				else
				{
					// fill Liquidity factor exposure
					double tmp_stom = isnormal(stom[date][symbol]) ? stom[date][symbol] : 0.;
					double tmp_stoq = isnormal(stoq[date][symbol]) ? stoq[date][symbol] : 0.;
					double tmp_stoa = isnormal(stoa[date][symbol]) ? stoa[date][symbol] : 0.;
					Liquidity[date][symbol] = 0.35 * tmp_stom + 0.35 * tmp_stoq + 0.3 * tmp_stoa;
				}
			}
		} // loop of symbols	

		// linear fit for nl-size
		if (vec_size.size() == vec_cubsize.size())
		{
			//double c0, c1, cov00, cov01, cov11, sumsq;
			//gsl_fit_linear(&vec_size[0], 1, &vec_cubsize[0], 1, vec_size.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
			double c1, cov11, sumsq;
			gsl_fit_mul(&vec_size[0], 1, &vec_cubsize[0], 1, vec_size.size(), &c1, &cov11, &sumsq);
			for (int s = 0; s != vec_size.size(); ++s)
			{
				//double res = vec_cubsize[s] - (c0 + c1 * vec_size[s]);
				double res = vec_cubsize[s] - (c1 * vec_size[s]);
				// fill NlSize factor exposure
				NlSize[date][vec_tmp[s]] = res;
			}			
		}

		// normalize the factor exposures
		// size, beta and momentum already normalize as descriptors
		// resvol and liquidity are composed by descriptors, need to re-normalize
		//normalize(Size[date], MarketData::instance().market_data[date].cap_weight);
		//normalize(Beta[date], MarketData::instance().market_data[date].cap_weight);
		//normalize(Momentum[date], MarketData::instance().market_data[date].cap_weight);
		normalize(ResVol[date], MarketData::instance().market_data[date].cap_weight);
		normalize(NlSize[date], MarketData::instance().market_data[date].cap_weight);
		normalize(Liquidity[date], MarketData::instance().market_data[date].cap_weight);

		// fillna for the exposures
		fillna_mean(Size[date], MarketData::instance().market_data[date].cap_weight);
		fillna_mean(Beta[date], MarketData::instance().market_data[date].cap_weight);
		fillna_mean(Momentum[date], MarketData::instance().market_data[date].cap_weight);
		fillna_mean(ResVol[date], MarketData::instance().market_data[date].cap_weight);
		fillna_mean(NlSize[date], MarketData::instance().market_data[date].cap_weight);
		fillna_mean(Liquidity[date], MarketData::instance().market_data[date].cap_weight);

		// check exposure values and retrieve value for ResVol fit
		vector<double> vec_resvol, vec_beta;
        vec_resvol.clear(); vec_beta.clear(); vec_size.clear(); vec_tmp.clear();
		for (int s = 0; s != MarketData::instance().symbols.size(); ++s)
		{
			string symbol = MarketData::instance().symbols[s];
			/*
			double tmp_size = Size[date][symbol];
			double tmp_beta = Beta[date][symbol];
			double tmp_resvol = ResVol[date][symbol];
			double tmp_momentum = Momentum[date][symbol];
			double tmp_nlsize = NlSize[date][symbol];
			double tmp_liquidity = Liquidity[date][symbol];
			if (!isnormal(tmp_size))
				cout<<"Abnormal Size exposure: "<<tmp_size<<endl;
			if (!isnormal(tmp_beta))
				cout<<"Abnormal Beta exposure: "<<tmp_beta<<endl;
			if (!isnormal(tmp_momentum))
				cout<<"Abnormal Momentum exposure: "<<tmp_momentum<<endl;
			if (!isnormal(tmp_resvol))
				cout<<"Abnormal ResVol exposure: "<<tmp_resvol<<endl;
			if (!isnormal(tmp_nlsize))
				cout<<"Abnormal NlSize exposure: "<<tmp_nlsize<<endl;
			if (!isnormal(tmp_liquidity))
				cout<<"Abnormal Liquidity exposure: "<<tmp_liquidity<<endl;
				*/

			vec_resvol.push_back(ResVol[date][symbol]);
			vec_beta.push_back(Beta[date][symbol]);
			vec_size.push_back(Size[date][symbol]);
			vec_tmp.push_back(symbol);
		}
	
		// remove collinearity between ResVol adn Size,Beta
		if (vec_size.size() == vec_beta.size() && vec_beta.size() == vec_resvol.size() && vec_size.size() != 0)
		{
			int n = vec_size.size();
			double chisq;
			gsl_matrix* X = gsl_matrix_alloc(n, 2);
			gsl_vector* Y = gsl_vector_alloc(n);
			gsl_vector* c = gsl_vector_alloc(2);
			gsl_matrix* cov = gsl_matrix_alloc(2, 2);
			for (int i = 0; i != n; ++i)
			{
				gsl_matrix_set(X, i, 0, vec_size[i]);
				gsl_matrix_set(X, i, 1, vec_beta[i]);
				gsl_vector_set(Y, i, vec_resvol[i]);				
			}

			gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(n, 2);
			gsl_multifit_linear(X, Y, c, cov, &chisq, work);
			gsl_multifit_linear_free(work);
			
			double c0 = gsl_vector_get(c, 0);
			double c1 = gsl_vector_get(c, 1);
			for (int i = 0; i != n; ++i)
				ResVol[date][vec_tmp[i]] = ((vec_resvol[i] - (c0 * vec_size[i] + c1 * vec_beta[i])));

			gsl_matrix_free(X);
  			gsl_vector_free(Y);
		  	gsl_vector_free(c);
		  	gsl_matrix_free(cov);
		}

		// normalize ResVol exposure again
		normalize(ResVol[date], MarketData::instance().market_data[date].cap_weight);
		fillna_mean(ResVol[date], MarketData::instance().market_data[date].cap_weight);
	} // loop of trading days

	auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Time taken by compute_exposure: "<<(duration.count() / 1000000.)<<" seconds"<<endl;
}

void StyleFactor::compute_returns()
{
	auto start = chrono::high_resolution_clock::now();
	// apply linear regression per day with the six style factors and the country factor
	for (int i = (504 + 21); (i + 10) != MarketData::instance().trading_days.size(); ++i)
	{
		int date = MarketData::instance().trading_days[i];
		int f_date = MarketData::instance().trading_days[i + 1];
		// retrieve the values
		vector<double> vec_rt, vec_size, vec_beta, vec_momentum, vec_resvol, vec_nlsize, vec_liquidity, vec_w;
		vec_rt.clear(); vec_size.clear(); vec_beta.size(); vec_momentum.clear(); vec_resvol.clear(); vec_nlsize.clear(); vec_liquidity.clear(); vec_w.clear();
		for (int s = 0; s != MarketData::instance().symbols.size(); ++s)
		{
			string symbol = MarketData::instance().symbols[s];
			// use return of next tradig day
			double rt = isnormal(MarketData::instance().market_data[f_date].returns[s]) ? MarketData::instance().market_data[f_date].returns[s] : MarketData::instance().market_data[f_date].mkt_return;
			double rft = MarketData::instance().market_data[f_date].shibor;
			double cap_w = isnormal(MarketData::instance().market_data[date].tot_mkt_val[s]) ? sqrt(MarketData::instance().market_data[date].tot_mkt_val[s]) : 0.;
			double tmp_size = Size[date][symbol];
			double tmp_beta = Beta[date][symbol];
			double tmp_momentum = Momentum[date][symbol];
			double tmp_resvol = ResVol[date][symbol];
			double tmp_nlsize = NlSize[date][symbol];
			double tmp_liquidity = Liquidity[date][symbol];
			
			if (!isnormal(rt)) cout<<"Abnormal rt value: "<<rt<<endl;

			vec_rt.push_back(rt - rft);
			vec_w.push_back(cap_w);
			vec_size.push_back(tmp_size);
			vec_beta.push_back(tmp_beta);
			vec_momentum.push_back(tmp_momentum);
			vec_resvol.push_back(tmp_resvol);
			vec_nlsize.push_back(tmp_nlsize);
			vec_liquidity.push_back(tmp_liquidity);
		} // loop of symbols

		// multi-var linear fit
		if (vec_rt.size() != 0)
		{
			int n = vec_rt.size();
            double chisq;
            gsl_matrix* X = gsl_matrix_alloc(n, 7);
            gsl_vector* Y = gsl_vector_alloc(n);
			gsl_vector* W = gsl_vector_alloc(n);
            gsl_vector* c = gsl_vector_alloc(7);
            gsl_matrix* cov = gsl_matrix_alloc(7, 7);

			for (int s = 0; s != n; ++s)
			{
				gsl_matrix_set(X, s, 0, 1.);
                gsl_matrix_set(X, s, 1, vec_size[s]);
				gsl_matrix_set(X, s, 2, vec_beta[s]);
				gsl_matrix_set(X, s, 3, vec_momentum[s]);
				gsl_matrix_set(X, s, 4, vec_resvol[s]);
				gsl_matrix_set(X, s, 5, vec_nlsize[s]);
				gsl_matrix_set(X, s, 6, vec_liquidity[s]);
                gsl_vector_set(Y, s, vec_rt[s]);
				gsl_vector_set(W, s, vec_w[s]);
			}

			gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(n, 7);
            gsl_multifit_wlinear(X, W, Y, c, cov, &chisq, work);
            gsl_multifit_linear_free(work);

			// factor returns
            double c1 = gsl_vector_get(c, 1);
			double c2 = gsl_vector_get(c, 2);
            double c3 = gsl_vector_get(c, 3);
			double c4 = gsl_vector_get(c, 4);
            double c5 = gsl_vector_get(c, 5);
			double c6 = gsl_vector_get(c, 6);

			// fill the values
			// skip bad fit
			if (!(fabs(c1) > 1. || fabs(c2) > 1. || fabs(c3) > 1. || fabs(c4) > 1. || fabs(c5) > 1. || fabs(c6) > 1.))
			{
				dates.push_back(date);
				Size_Rt[date] = c1;
				Beta_Rt[date] = c2;
				Momentum_Rt[date] = c3;
				ResVol_Rt[date] = c4;
				NlSize_Rt[date] = c5;
				Liquidity_Rt[date] = c6;
			}

			gsl_matrix_free(X);
            gsl_vector_free(Y);
			gsl_vector_free(W);
            gsl_vector_free(c);
            gsl_matrix_free(cov);
		} // WLS fit
	} // loop of trading days

	// save to vectors
	for (int i = 0; i != dates.size(); ++i)
	{
		int date = dates[i];
		vector<double> tmp_vec;
		tmp_vec.clear();
		// push back factor returns
		tmp_vec.push_back(Size_Rt[date]);
		tmp_vec.push_back(Beta_Rt[date]);
		tmp_vec.push_back(Momentum_Rt[date]);
		tmp_vec.push_back(ResVol_Rt[date]);
		tmp_vec.push_back(NlSize_Rt[date]);
		tmp_vec.push_back(Liquidity_Rt[date]);
		
		factor_rt.push_back(tmp_vec);
	}

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Time taken by compute_returns: "<<(duration.count() / 1000000.)<<" seconds"<<endl;
}
