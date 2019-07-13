// Created by Ruiqi Zhang on 2019-07-08

#ifndef MARKETDATA_H
#define MARKETDATA_H

#include "singleton.h"
#include <vector>
#include <map>

// market data structure
struct Data
{
	int date;
	double shibor;
	std::vector<double> shares;
	std::vector<double> volume;
	std::vector<double> tot_mkt_val;
	std::vector<double> cap_weight;
	std::vector<double> turnover;
	std::vector<double> returns;
	double mkt_return;
};

// singleton mode allowing for only one instance
class MarketData : public Singleton<MarketData>
{
	public:
		MarketData() = default;
		~MarketData() final = default;

		void init_market_data();

	public:
		std::vector<int> trading_days;
		std::vector<std::string> symbols;
		std::map<int, Data> market_data;

		inline double safe_stod(std::string str);
};

#endif // MARKETDATA_H
