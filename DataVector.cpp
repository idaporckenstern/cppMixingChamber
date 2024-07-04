#include "DataVector.h"
#include <iostream>

DataVector::DataVector(std::vector<double> xVector, std::vector<double> yVector)
{
	this->xVector = xVector;
	this->yVector = yVector;
}

void DataVector::initDataVector()
{
	//initializes the data vector filled with zeros so they can be changed later

	this->data.resize(this->xVector.size());
	for (int i = 0; i < this->xVector.size(); ++i)
	{
		this->data[i].resize(this->yVector.size());
	}

}

void DataVector::setData(double dataPoint, int i, int j)
{
	this->data[i][j] = dataPoint;
}

double DataVector::getData(int i, int j)
{
	return this->data[i][j];
}

double DataVector::getXPoint(int i)
{
	return this->xVector[i];
}

double DataVector::getYPoint(int j)
{
	return this->yVector[j];
}

void DataVector::testing()
{
	for (int i = 0; i < data.size(); ++i)
	{
		for (int j = 0; j < data[i].size(); ++j)
		{
			std::cout << data[i][j];
		}
		std::cout << "\n";
	}
}