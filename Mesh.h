#pragma once
#include <vector>
#include "DataVector.h"
#include <memory>
#include "Opening.h"
#include "Outlet.h"
#include <iostream>
#include <cmath>

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
		//the first two i and j loops set the boundary conditions for u
		for (int i = 0; i < this->M + 1; ++i)
		{
			this->uVelocity->setData(-uVelocity->getData(i, 1), i, 0);
			this->uVelocity->setData(-uVelocity->getData(i, this->N), i, this->N + 1);

			for (Outlet outlet : outlets)
			{
				if ((this->uVelocity->getXPoint(i) >= outlet.getStartingPoint().x) && (this->uVelocity->getXPoint(i) <= outlet.getEndPoint().x) && outlet.getIsOpen())
				{
					if (outlet.getStartingPoint().y == 0)
					{
						this->uVelocity->setData(this->uVelocity->getData(i, 1), i, 0);
					}
					else
					{
						this->uVelocity->setData(this->uVelocity->getData(i, N + 1), i, N);
					}
				}
			}
		}

		for (int j = 0; j < this->N + 2; ++j)
		{
			this->uVelocity->setData(0, 0, j);
			this->uVelocity->setData(0, M, j);

			for (int k = 0; k < inletConstants.size(); ++k)
			{
				if ((this->uVelocity->getYPoint(j) >= inlets[k].getStartingPoint().y) && (this->uVelocity->getYPoint(j) <= inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == 0))
				{
					double z = uVelocity->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->uVelocity->setData(inletConstants[k] * lambda(z, this->lengthY / 8), 0, j);
				}
				if ((this->uVelocity->getYPoint(j) >= inlets[k].getStartingPoint().y) && (uVelocity->getYPoint(j) <= inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == this->lengthX))
				{
					double z = this->uVelocity->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->uVelocity->setData(-inletConstants[k] * lambda(z, this->lengthY / 8), M, j);
				}
			}
		}

		//the next two sets of loops set up the boundary conditions for v
		for (int i = 0; i < this->M + 2; ++i)
		{
			this->vVelocity->setData(0, i, 0);
			this->vVelocity->setData(0, i, N);

			for (Outlet outlet : this->outlets)
			{
				if ((this->vVelocity->getXPoint(i) >= outlet.getStartingPoint().x) && (this->vVelocity->getXPoint(i) <= outlet.getEndPoint().x) && outlet.getIsOpen())
				{
					if (outlet.getStartingPoint().y == 0)
					{
						this->vVelocity->setData((4 * this->vVelocity->getData(i, 1) - this->vVelocity->getData(i, 2)) / 3, i, 0);
					}
					else
					{
						this->vVelocity->setData((4 * this->vVelocity->getData(i, this->N - 1) - this->vVelocity->getData(i, this->N - 2)) / 3, i, N);
					}
				}
			}

			for (int k = 0; k < inletConstants.size(); ++k)
			{
				if ((this->vVelocity->getXPoint(i) >= inlets[k].getStartingPoint().x) && (this->vVelocity->getXPoint(i) <= inlets[k].getEndPoint().x) && (inlets[k].getStartingPoint().y == 0))
				{
					double z = vVelocity->getXPoint(i) - inlets[k].getStartingPoint().x;
					this->vVelocity->setData(inletConstants[k] * lambda(z, this->lengthX / 8), i, 0);
				}
				if ((this->vVelocity->getXPoint(i) > inlets[k].getStartingPoint().x) && (vVelocity->getXPoint(i) < inlets[k].getEndPoint().x) && (inlets[k].getStartingPoint().y == this->lengthY))
				{
					double z = this->vVelocity->getXPoint(i) - inlets[k].getStartingPoint().x;
					this->vVelocity->setData(-inletConstants[k] * lambda(z, this->lengthX / 8), i, N);
				}
			}
		}

		for (int j = 0; j < this->N + 1; ++j)
		{
			this->vVelocity->setData(-this->vVelocity->getData(1, j), 0, j);
			this->vVelocity->setData(-this->vVelocity->getData(M, j), M + 1, j);

			for (int k = 0; k < inletConstants.size(); ++k)
			{
				
				if ((this->vVelocity->getYPoint(j) >= inlets[k].getStartingPoint().y) && (this->vVelocity->getYPoint(j) <= inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == 0))
				{
					double z = vVelocity->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->vVelocity->setData(2 * (inletConstants[k] * lambda(z, this->lengthY / 8)) - this->vVelocity->getData(1, j), 0, j);
				}
				if ((this->vVelocity->getYPoint(j) >= inlets[k].getStartingPoint().y) && (vVelocity->getYPoint(j) <= inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == this->lengthX))
				{
					double z = this->vVelocity->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->vVelocity->setData(2 * (-inletConstants[k] * lambda(z, this->lengthY / 8)) - this->vVelocity->getData(M, j), M + 1, j);
				}
			}
		}

		//the next two sets of loops set up the boundary conditions for Y
		int fY;
		if (std::fmod(this->t, 10) < 5)
		{
			fY = 1;
		}
		else
		{
			fY = 0;
		}
		
		for (int i = 0; i < M + 2; ++i)
		{
			this->Y->setData(this->Y->getData(i, 1), i, 0);
			this->Y->setData(this->Y->getData(i, N), i, N + 1);

			for (int k = 0; k < inletConstants.size(); ++k)
			{
				if ((this->Y->getXPoint(i) >= inlets[k].getStartingPoint().x) && (this->Y->getXPoint(i) <= inlets[k].getEndPoint().x) && (inlets[k].getStartingPoint().y == 0))
				{
					double z = Y->getXPoint(i) - inlets[k].getStartingPoint().x;
					this->Y->setData(2 * fY - this->Y->getData(i, 1), i, 0);
					
				}
				if ((this->Y->getXPoint(i) > inlets[k].getStartingPoint().x) && (Y->getXPoint(i) < inlets[k].getEndPoint().x) && (inlets[k].getStartingPoint().y == this->lengthY))
				{
					double z = this->Y->getXPoint(i) - inlets[k].getStartingPoint().x;
					this->Y->setData(2 * fY - this->Y->getData(i, N), i, N + 1);
				}
			}
		}

		for (int j = 0; j < N + 2; ++j)
		{
			this->Y->setData(this->Y->getData(1, j), 0, j);
			this->Y->setData(this->Y->getData(M, j), M + 1, j);

			for (int k = 0; k < inletConstants.size(); ++k)
			{

				if ((this->Y->getYPoint(j) >= inlets[k].getStartingPoint().y) && (this->Y->getYPoint(j) <= inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == 0))
				{
					double z = Y->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->Y->setData(2 * fY - this->Y->getData(1, j), 0, j);
				}
				if ((this->Y->getYPoint(j) >= inlets[k].getStartingPoint().y) && (Y->getYPoint(j) <= inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == this->lengthX))
				{
					double z = this->Y->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->Y->setData(2 * fY - this->Y->getData(M, j), M + 1, j);
				}
			}
		}
	}
};
