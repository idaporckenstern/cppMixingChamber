#pragma once
#include <vector>
#include "DataVector.h"
#include <memory>
#include "Opening.h"
#include "Outlet.h"
#include <iostream>

class Mesh
{
private:
	//initial variables
	double t;
	double tEnd;
	double M; //amount of points in the mesh on the x direction
	double N; //amount of points in the mesh on the y direction
	double lengthX;
	double lengthY;
	double h;
	double Re;
	double Sc;

	//inlets and outlets
	std::vector<Opening> inlets;
	std::vector<Outlet> outlets;


	//wanted data
	std::unique_ptr<DataVector> uVelocity;
	std::unique_ptr<DataVector> vVelocity;
	std::unique_ptr<DataVector> Y;
	std::unique_ptr<DataVector> pressure;

	//double data for FTCS
	std::unique_ptr<std::vector<double>> u2;
	std::unique_ptr<std::vector<double>> v2;
	std::unique_ptr<std::vector<double>> Y2;

public:
	Mesh(double M, double N, double lengthX, double lengthY, double Re, double Sc, double t, double tEnd, std::vector<Opening>& inlets, std::vector<Outlet>& outlets);
	void testing();

	double getT() { return this->t; }
	double getTEnd() { return this->tEnd; }

	template<typename T>
	void setBoundaryConditions(T lambda, std::vector<double>& inletConstants)
	{
		for (int i = 0; i < this->M + 1; ++i)
		{
			this->uVelocity->setData(-uVelocity->getData(i, 1), i, 0);
			this->uVelocity->setData(-uVelocity->getData(i, this->N), i, this->N + 1);

			for (Outlet outlet : outlets)
			{
				if (this->uVelocity->getXPoint(i) >= outlet.getStartingPoint().x)
				{

				}
			}
		}
	}
};
