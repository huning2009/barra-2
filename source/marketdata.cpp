// Created by Ruiqi Zhang on 2019-07-08

#include "marketdata.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include<boost/algorithm/string.hpp>
using namespace std;

double MarketData::safe_stod(string str)
{
	try
	{
		return stod(str);
	}
	catch(...)
	{
		return numeric_limits<double>::quiet_NaN();
	}
}

void MarketData::init_market_data()
{
	auto start = chrono::high_resolution_clock::now();
	// open input files
	ifstream f_price("/Users/ruiqizhang/workarea/barra/data/AdjClose.csv");
	ifstream f_shares("/Users/ruiqizhang/workarea/barra/data/Shares.csv");
	ifstream f_shibor("/Users/ruiqizhang/workarea/barra/data/Shibor3M.csv");
	ifstream f_tot_mkt_val("/Users/ruiqizhang/workarea/barra/data/TotalMarketValue.csv");
	ifstream f_vol("/Users/ruiqizhang/workarea/barra/data/Volume.csv");
	ifstream f_turnover("/Users/ruiqizhang/workarea/barra/data/TurnoverRate.csv");
	ifstream f_trade_stk("/Users/ruiqizhang/workarea/barra/data/trade_stocks.csv");

	int N = 1500; // maximum dates
	int S = 1000; // maximum symbols
	// lines to record content
	string s_price, s_price_prev, s_shares, s_shibor, s_tot_mkt_val, s_vol, s_turnover, s_trade_stk;

	// 2841 lines in these files
	vector<double> vec_price, vec_price_prev;
	vec_price.clear(); vec_price_prev.clear();
	int l = 0;
	while (l < N)
	//while (l < 1000)
	{
		// market data of this date
		Data data;

		getline(f_price, s_price);
		getline(f_shares, s_shares);
		getline(f_shibor, s_shibor);
		getline(f_tot_mkt_val, s_tot_mkt_val);
		getline(f_vol, s_vol);
		getline(f_turnover, s_turnover);
		getline(f_trade_stk, s_trade_stk);

		// init symbols using the first line of trade_stocks
		if (l == 0)
		{
			vector<string> out_symbols;
			out_symbols.clear();
			boost::split(out_symbols, s_trade_stk, boost::is_any_of(","));
			// need to drop first empty elemnt 
			for (int i = 1; i != S; ++i)
			//for (int i = 1; i != 500; ++i)
				symbols.push_back(out_symbols[i]);
			// then skip the first line
			l++;
			continue;
		}

		vector<string> out_price, out_shares, out_shibor, out_tot_mkt_val, out_vol, out_turnover, out_trade_stk;
		out_price.clear(); out_shares.clear(); out_shibor.clear(); out_tot_mkt_val.clear(); out_vol.clear(); out_turnover.clear(); out_trade_stk.clear();
		boost::split(out_price, s_price, boost::is_any_of(","));
		boost::split(out_shares, s_shares, boost::is_any_of(","));
		boost::split(out_shibor, s_shibor, boost::is_any_of(","));
		boost::split(out_tot_mkt_val, s_tot_mkt_val, boost::is_any_of(","));
		boost::split(out_vol, s_vol, boost::is_any_of(","));
		boost::split(out_turnover, s_turnover, boost::is_any_of(","));
		boost::split(out_trade_stk, s_trade_stk, boost::is_any_of(","));

		/*
		// check size of the lines
		if (out_price.size() != 3684) { cout<<"out_price at "<<l<<" line, "<<out_price.size()<<endl; }
		if (out_shares.size() != 3684) { cout<<"out_shares at "<<l<<" line, "<<out_shares.size()<<endl; }
		if (out_tot_mkt_val.size() != 3684) { cout<<"out_tot_mkt_val at "<<l<<" line, "<<out_tot_mkt_val.size()<<endl; }
		if (out_vol.size() != 3684) { cout<<"out_vol at "<<l<<" line, "<<out_vol.size()<<endl; }
		if (out_trade_stk.size() != 3684) { cout<<"out_trade_stk at "<<l<<" line, "<<out_trade_stk.size()<<endl; }
		*/

		// skip first day due to return calculation
		// retrive valid trading dates from Shibor3M.csv
		if (l != 1)
		{
			trading_days.push_back(stoi(out_shibor[0]));
			data.date = stoi(out_shibor[0]);
			data.shibor = 0.01 * safe_stod(out_shibor[1]) / 63.;
		}

		// 3683 is the length of the tables. i.e. num of stocks
		vector<double> vec_tmp_price;
		vec_tmp_price.clear();
		// cap weight
		vector<double> vec_cap;
		vec_cap.clear();
		double tot_cap = 0.;
		for (int i = 1; i != S; ++i)
		//for (int i = 1; i != 500; ++i)
		{
			// note NaNs
			/*
			double tmp_shares = isnormal(out_shares[i]) ? safe_stod(out_shares[i]) : numeric_limits<double>::quiet_NaN();
			double tmp_tot_mkt_val = isnormal(out_tot_mkt_val[i]) ? safe_stod(out_tot_mkt_val[i]) : numeric_limits<double>::quiet_NaN();
			double tmp_volume = isnormal(out_vol[i]) ? safe_stod(out_vol[i]) : numeric_limits<double>::quiet_NaN();
			double tmp_price = isnormal(out_price[i]) ? safe_stod(out_price[i]) : numeric_limits<double>::quiet_NaN();
			*/
			double tmp_shares = safe_stod(out_shares[i]);
			double tmp_tot_mkt_val = safe_stod(out_tot_mkt_val[i]);
			double tmp_volume = safe_stod(out_vol[i]);
			double tmp_turnover = safe_stod(out_turnover[i]);
			double tmp_price = safe_stod(out_price[i]);
		
			if (isnormal(tmp_tot_mkt_val))
			{
				vec_cap.push_back(tmp_tot_mkt_val);
				tot_cap += tmp_tot_mkt_val;
			}
			else
				vec_cap.push_back(0.);

			if (l != 1)
			{
				data.shares.push_back(tmp_shares);
				data.tot_mkt_val.push_back(tmp_tot_mkt_val);
				data.volume.push_back(tmp_volume);
				data.turnover.push_back(tmp_turnover);
			}
			vec_tmp_price.push_back(tmp_price);
		}

		//cap weight
		if (tot_cap > 0.00001)
		{
			for (int s = 0; s != vec_cap.size(); ++s) 
				vec_cap[s] = vec_cap[s] / tot_cap;
		}
		data.cap_weight = vec_cap;

		// compute daily return
		if (l == 1)
			vec_price = vec_tmp_price;
		else
		{
			vec_price_prev = vec_price;
			vec_price = vec_tmp_price;
			for (int i = 1; i != vec_price.size(); ++i)
			{
				double tmp_rt = (isnormal(vec_price[i]) && isnormal(vec_price_prev[i])) ? ((vec_price[i] / vec_price_prev[i]) - 1.) : numeric_limits<double>::quiet_NaN();
				data.returns.push_back(tmp_rt);
			}
		}

		if (l == 1) { l++; continue; }
		// compute the market return
		double tmp_mkt_rt = 0.;
		double market_val_sum = 0;
		for (int i = 0; i != data.returns.size(); ++i)
			if (isnormal(data.tot_mkt_val[i]) && isnormal(data.returns[i]))
				market_val_sum += data.tot_mkt_val[i];
		for (int i = 0; i != data.returns.size(); ++i)
			tmp_mkt_rt += (isnormal(data.tot_mkt_val[i]) && isnormal(data.returns[i])) ? (data.tot_mkt_val[i] / market_val_sum) * data.returns[i] : 0.;
		data.mkt_return = tmp_mkt_rt;

		market_data[data.date] = data;	

		l++;
	} // loop of trading days

	f_price.close();
	f_shares.close();
	f_shibor.close();
	f_tot_mkt_val.close();
	f_vol.close();
	f_turnover.close();
	f_trade_stk.close();

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	cout << "Time taken by init_market_data: "<<(duration.count() / 1000000.)<<" seconds"<<endl;
}

