#pragma once
#include <vector>
#include <cmath>
#include "Functions.h"

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
	void replaceDataVector(std::vector<std::vector<double>> data);
	double getData(int i, int j);
	double getXPoint(int i);
	double getYPoint(int j);
	std::vector<std::vector<double>> getData() { return data; }
	double absoluteMax();
	double sumX(int column, int start, int end);
	double sumY(int row, int start, int end);
	void testing();
};
