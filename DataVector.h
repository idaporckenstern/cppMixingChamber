#pragma once
#include <vector>

class DataVector
{
private:
	std::vector<double> xVector;
	std::vector<double> yVector;

	std::vector<std::vector<double>> data;

public:
	DataVector(std::vector<double> xVector, std::vector<double> yVector);
	void initDataVector();

	void setData(double dataPoint, int i, int j);
	double getData(int i, int j);
	double getXPoint(int i);
	double getYPoint(int j);
	void testing();
};
