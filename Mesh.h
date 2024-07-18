#pragma once
#include <vector>
#include "DataVector.h"
#include <memory>
#include "Opening.h"
#include "Outlet.h"
#include <iostream>
#include <cmath>
#include<iterator>

class Mesh
{
private:
	//initial variables
	double t;
	double dt;
	double tEnd;
	double M; //amount of points in the mesh on the x direction
	double N; //amount of points in the mesh on the y direction
	double lengthX;
	double lengthY;
	double h;
	double Re;
	double Sc;
	bool doPlot = false;

	//inlets and outlets
	std::vector<Opening> inlets;
	std::vector<Outlet> outlets;


	//wanted data
	std::unique_ptr<DataVector> uVelocity;
	std::unique_ptr<DataVector> vVelocity;
	std::unique_ptr<DataVector> Y;
	std::unique_ptr<DataVector> pressure;
	std::unique_ptr<DataVector> R;

	//double data for FTCS
	std::unique_ptr<DataVector> u2;
	std::unique_ptr<DataVector> v2;
	std::unique_ptr<DataVector> Y2;

	//vectors to find dt
	std::vector<std::vector<double>> uHat;
	std::vector<std::vector<double>> vHat;

	//right hand side of the pressure equation
	std::unique_ptr<DataVector> pressureRightHandSide;

public:
	Mesh(double M, double N, double lengthX, double lengthY, double Re, double Sc, double t, double tEnd, std::vector<Opening>& inlets, std::vector<Outlet>& outlets);
	void testing();

	double getT() { return this->t; }
	double getTEnd() { return this->tEnd; }
	double getM() { return this->M; }
	double getN() { return this->N; }
	double getH() { return this->h; }
	double getUData(int i, int j) { return this->uVelocity->getData(i, j); }
	double getVData(int i, int j) { return this->vVelocity->getData(i, j); }
	double getDT() { return this->dt; }
	void setT(double t) { this->t = t; }
	void setDoPlot(bool doPlot) { this->doPlot = doPlot; }

	double solveDT(double CFL);
	double absoluteMax(std::vector<std::vector<double>> &vector);
	double min(double a, double b);
	

	void stepForward();
	void correctVelocities();
	void conservationMassCorrection();
	void solvePressure();
	std::vector<std::vector<double>> poissonSolver(int maxIterations, double epsilon);
	std::vector<std::vector<double>> multiGridSolver(std::vector<std::vector<double>> &phi, std::vector<std::vector<double>> &rightHandSide, double h);
	std::vector<std::vector<double>> centerGS(std::vector<std::vector<double>> &phi, std::vector<std::vector<double>> &rightHandSide, double h, int numIterations);
	std::vector<std::vector<double>> centerResidual(std::vector<std::vector<double>> &phi, std::vector<std::vector<double>> &rightHandSide, double h);
	std::vector<std::vector<double>> centerRestrictCells(std::vector<std::vector<double>> &vector);
	std::vector<std::vector<double>> centerProlongCells(std::vector<std::vector<double>> &vector);
	std::vector<std::vector<double>> addVectors(std::vector<std::vector<double>> &vector1, std::vector<std::vector<double>> &vector2);
	void setGSBoundaryConditions(std::vector<std::vector<double>> &phi);

	template<typename T>
	void setBoundaryConditionsU(T lambda, std::vector<double>& inletConstants, double t)
	{
		//the first two i and j loops set the boundary conditions for u
		for (int i = 0; i < this->M + 1; ++i)
		{
			this->uVelocity->setData(-uVelocity->getData(i, 0), i, 1);
			this->uVelocity->setData(-uVelocity->getData(i, this->N + 1), i, this->N);

			for (int k = 0; k < outlets.size(); ++k)
			{
				if ((this->uVelocity->getXPoint(i) >= outlets[k].getStartingPoint().x) && (this->uVelocity->getXPoint(i) <= outlets[k].getEndPoint().x) && outlets[k].getIsOpen())
				{
					if (outlets[k].getStartingPoint().y == 0)
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
				if ((this->uVelocity->getYPoint(j) > inlets[k].getStartingPoint().y) && (this->uVelocity->getYPoint(j) < inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == 0))
				{
					double z = uVelocity->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->uVelocity->setData(inletConstants[k] * lambda(z, this->lengthY / 8), 0, j);
				}
				if ((this->uVelocity->getYPoint(j) > inlets[k].getStartingPoint().y) && (uVelocity->getYPoint(j) < inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == this->lengthX))
				{
					double z = this->uVelocity->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->uVelocity->setData(-inletConstants[k] * lambda(z, this->lengthY / 8), M, j);
				}
			}
		}
	}

	template<typename T>
	void setBoundaryConditionsV(T lambda, std::vector<double>& inletConstants, double t)
	{

		//the next two sets of loops set up the boundary conditions for v
		for (int i = 0; i < this->M + 2; ++i)
		{
			this->vVelocity->setData(0, i, 0);
			this->vVelocity->setData(0, i, N);

			for (int k = 0; k < outlets.size(); ++k)
			{
				if ((this->vVelocity->getXPoint(i) > outlets[k].getStartingPoint().x) && (this->vVelocity->getXPoint(i) < outlets[k].getEndPoint().x) && outlets[k].getIsOpen())
				{
					if (outlets[k].getStartingPoint().y == 0)
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
				if ((this->vVelocity->getXPoint(i) > inlets[k].getStartingPoint().x) && (this->vVelocity->getXPoint(i) < inlets[k].getEndPoint().x) && (inlets[k].getStartingPoint().y == 0))
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

				if ((this->vVelocity->getYPoint(j) > inlets[k].getStartingPoint().y) && (this->vVelocity->getYPoint(j) < inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == 0))
				{
					double z = vVelocity->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->vVelocity->setData(2 * (inletConstants[k] * lambda(z, this->lengthY / 8)) - this->vVelocity->getData(1, j), 0, j);
				}
				if ((this->vVelocity->getYPoint(j) > inlets[k].getStartingPoint().y) && (vVelocity->getYPoint(j) < inlets[k].getEndPoint().y) && (inlets[k].getStartingPoint().x == this->lengthX))
				{
					double z = this->vVelocity->getYPoint(j) - inlets[k].getStartingPoint().y;
					this->vVelocity->setData(2 * (-inletConstants[k] * lambda(z, this->lengthY / 8)) - this->vVelocity->getData(M, j), M + 1, j);
				}
			}
		}
	}

	template<typename T>
	void setBoundaryConditionsY(T lambda, std::vector<double>& inletConstants, double t)
	{
		//the next two sets of loops set up the boundary conditions for Y
		int fY;
		if (std::fmod(t, 10) < 5)
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
