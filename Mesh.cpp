#include "Mesh.h"

Mesh::Mesh(double M, double N, double lengthX, double lengthY, double Re, double Sc, double t, double tEnd, std::vector<Opening> &inlets, std::vector<Outlet> &outlets)
	: M(M), N(N), lengthX(lengthX), lengthY(lengthY), h(lengthX/M), Re(Re), Sc(Sc), t(t), tEnd(tEnd),
	uVelocity(std::make_unique<DataVector>(linspace(0.0, lengthX, M + 1), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	u2(std::make_unique<DataVector>(linspace(0.0, lengthX, M + 1), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	vVelocity(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0, lengthY, N + 1))),
	v2(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0, lengthY, N + 1))),
	Y(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	Y2(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	pressure(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	R(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	pressureRightHandSide(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	inlets(inlets), outlets(outlets)

{
	uVelocity->initDataVector();
	u2->initDataVector();
	vVelocity->initDataVector();
	v2->initDataVector();
	Y->initDataVector();
	Y2->initDataVector();
	pressure->initDataVector();
	R->initDataVector();
	pressureRightHandSide->initDataVector();

	this->uHat.resize(this->M);
	for (int i = 0; i < this->N; ++i)
	{
		this->uHat[i].resize(this->N);
	}

	this->vHat.resize(this->M);
	for (int i = 0; i < this->N; ++i)
	{
		this->vHat[i].resize(this->N);
	}
}

void Mesh::testing()
{
	
	//std::cout << "U2: \n";
	//u2->testing();
	std::cout << "U: \n";
	uVelocity->testing();
	std::cout << std::endl;
	std::cout << "V: \n";
	vVelocity->testing();
	//std::cout << "v2: \n";
	//v2->testing();
	std::cout << std::endl;
	//std::cout << "Y: \n";
	//Y->testing();
	//std::cout << std::endl;
	//std::cout << "Y2: \n";
	//Y2->testing();
	//std::cout << std::endl;
	std::cout << "pressure: \n";
	pressure->testing();
	std::cout << std::endl;
	//std::cout << "pressure rhs: \n";
	//pressureRightHandSide->testing();
	//std::cout << std::endl;
	//*/

	//std::cout << this->pressure->getData(1, 1);
}

double Mesh::solveDT(double CFL)
{
	double dtU = 0.0;
	double dtV = 0.0;
	double dtY = 0.0;
	double dtPara = 0.0;

	if (this->t == 0.0)
	{
		for (int j = 1; j < this->N + 1; ++j)
		{
			for (int i = 0; i < this->M; ++i)
			{
				this->uHat[i][j - 1] = (this->uVelocity->getData(i, j) + this->uVelocity->getData(i + 1, j)) / 2;
				
			}
		}

		for (int j = 0; j < this->N; ++j)
		{
			for (int i = 1; i < this->M + 1; ++i)
			{
				this->vHat[i - 1][j] = (this->vVelocity->getData(i, j) + this->vVelocity->getData(i, j + 1)) / 2;
			}
		}
	}

	dtU = CFL * this->h / (2 * this->uVelocity->absoluteMax() + this->vVelocity->absoluteMax());
	dtV = CFL * this->h / (2 * this->vVelocity->absoluteMax() + this->uVelocity->absoluteMax());
	dtY = CFL * this->h / (this->absoluteMax(this->uHat) + this->absoluteMax(this->vHat));
	dtPara = CFL * this->Re * std::pow(this->h, 2) / 4;
	this->dt = this->min(this->min(dtU, dtV), this->min(dtY, dtPara));
	return dt;
}

double Mesh::absoluteMax(std::vector<std::vector<double>> &vector)
{
	double max = 0.0;

	for (int i = 0; i < vector.size(); ++i)
	{
		for (int j = 0; j < vector[0].size(); ++j)
		{
			if (max < std::abs(vector[i][j]))
			{

				max = std::abs(vector[i][j]);
			}
		}
	}
	return max;
}

double Mesh::min(double a, double b)
{
	if (a < b)
	{
		return a;
	}
	return b;
}

void Mesh::stepForward()
{
	//solving for u^n+1
	for (int j = 1; j < this->N + 1; ++j)
	{
		for (int i = 1; i < this->M; ++i)
		{
			this->u2->setData(this->uVelocity->getData(i, j) + this->dt * (1 / this->Re * (((this->uVelocity->getData(i + 1, j) - 2 * this->uVelocity->getData(i, j) + this->uVelocity->getData(i - 1, j)) / (std::pow(this->h, 2))) + ((this->uVelocity->getData(i, j + 1) - 2 * this->uVelocity->getData(i, j) + this->uVelocity->getData(i, j - 1)) / (std::pow(this->h, 2)))) -
				((std::pow(((this->uVelocity->getData(i + 1, j) + this->uVelocity->getData(i, j)) / 2), 2) - std::pow(((this->uVelocity->getData(i, j) + this->uVelocity->getData(i - 1, j)) / 2), 2)) / this->h) -
				((((this->vVelocity->getData(i, j) + this->vVelocity->getData(i + 1, j)) / 2) * ((this->uVelocity->getData(i, j) + this->uVelocity->getData(i, j + 1)) / 2) - ((this->vVelocity->getData(i, j - 1) + this->vVelocity->getData(i + 1, j - 1)) / 2) * ((this->uVelocity->getData(i, j) + this->uVelocity->getData(i, j - 1)) / 2)) / this->h)), i, j);
		}
	}

	//solving for v^n + 1
	for (int j = 1; j < this->N; ++j)
	{
		for (int i = 1; i < this->M + 1; ++i)
		{
			this->v2->setData(this->vVelocity->getData(i, j) + this->dt * (1 / this->Re * (((this->vVelocity->getData(i + 1, j) - 2 * this->vVelocity->getData(i, j) + this->vVelocity->getData(i - 1, j)) / (std::pow(this->h, 2)))
				+ ((this->vVelocity->getData(i, j + 1) - 2 * this->vVelocity->getData(i, j) + this->vVelocity->getData(i, j - 1)) / (std::pow(this->h, 2)))) - ((((this->uVelocity->getData(i, j) + this->uVelocity->getData(i, j + 1)) / 2) * ((this->vVelocity->getData(i, j) + this->vVelocity->getData(i + 1, j)) / 2)
					- ((this->uVelocity->getData(i - 1, j) + this->uVelocity->getData(i - 1, j + 1)) / 2) * ((this->vVelocity->getData(i, j) + this->vVelocity->getData(i - 1, j)) / 2)) / this->h) - ((std::pow(((this->vVelocity->getData(i, j) + this->vVelocity->getData(i, j + 1)) / 2), 2)
						- std::pow(((this->vVelocity->getData(i, j) + this->vVelocity->getData(i, j - 1)) / 2), 2)) / this->h)), i, j);
		}
	}

	//solving for Y^n + 1
	for (int j = 1; j < this->N + 1; ++j)
	{
		for (int i = 1; i < this->M + 1; ++i)
		{
			this->Y2->setData(this->Y->getData(i, j) + this->dt * (1 / (this->Re * this->Sc) * (((this->Y->getData(i + 1, j) - 2 * this->Y->getData(i, j) + this->Y->getData(i - 1, j)) / (std::pow(this->h, 2)))
				+ ((this->Y->getData(i, j + 1) - 2 * this->Y->getData(i, j) + this->Y->getData(i, j - 1)) / (std::pow(this->h, 2))))) + (-((this->uVelocity->getData(i, j) + this->uVelocity->getData(i - 1, j)) / 2 * this->dt / this->h) / 2 * (this->Y->getData(i + 1, j)
					- this->Y->getData(i - 1, j)) + std::pow(((this->uVelocity->getData(i, j) + this->uVelocity->getData(i - 1, j)) / 2 * this->dt / this->h), 2) / 2 * (this->Y->getData(i + 1, j) - 2 * this->Y->getData(i, j) + this->Y->getData(i - 1, j))) + 
				(-((this->vVelocity->getData(i, j) + this->vVelocity->getData(i, j - 1)) / 2 * this->dt / this->h) / 2 * (this->Y->getData(i, j + 1) - this->Y->getData(i, j - 1)) + std::pow(((this->vVelocity->getData(i, j) + this->vVelocity->getData(i, j - 1)) / 2 * this->dt / this->h), 2) / 2 * (this->Y->getData(i, j + 1) 
					- 2 * this->Y->getData(i, j) + this->Y->getData(i, j - 1))), i, j);
		}
	}

	//move new data onto off of our copies and apply boundary conditions
	this->uVelocity->replaceDataVector(this->u2->getData());
	this->vVelocity->replaceDataVector(this->v2->getData());
	this->Y->replaceDataVector(this->Y2->getData());
}

void Mesh::correctVelocities()
{
	conservationMassCorrection();
	solvePressure();
}

void Mesh::conservationMassCorrection()
{
	double volumeFlow = this->uVelocity->sumY(0, 1, this->N + 1) * this->h - this->uVelocity->sumY(this->M, 1, this->N + 1) * this->h + this->vVelocity->sumX(0, 1, this->M + 1) * this->h - this->vVelocity->sumX(this->N, 1, this->M + 1) * this->h;

	double outletsOpen = 0;
	double openRatio;
	for (int i = 0; i < this->outlets.size(); ++i)
	{
		if (outlets[i].getIsOpen())
		{
			++outletsOpen;
		}	
	}
	openRatio = outletsOpen / outlets.size();
	double velocityCorrection = volumeFlow / openRatio;

	for (int i = 0; i < this->M + 1; ++i)
	{
		for (int k = 0; k < outlets.size(); ++k)
		{
			if ((this->vVelocity->getXPoint(i) >= outlets[k].getStartingPoint().x) && (this->vVelocity->getXPoint(i) <= outlets[k].getEndPoint().x) && outlets[k].getIsOpen())
			{
				if (outlets[k].getStartingPoint().y == 0)
				{
					this->vVelocity->setData(this->vVelocity->getData(i, 0) - velocityCorrection, i, 0);
				}
				else
				{
					this->vVelocity->setData(this->vVelocity->getData(i, this->N) + velocityCorrection , i, N);
				}
			}
		}
	}
}

void Mesh::solvePressure()
{
	for (int j = 1; j < this->N + 1; ++j)
	{
		for (int i = 1; i < this->M + 1; ++i)
		{
			this->pressureRightHandSide->setData(1 / this->dt * ((this->uVelocity->getData(i, j) - this->uVelocity->getData(i - 1, j)) / this->h + (this->vVelocity->getData(i, j) - this->vVelocity->getData(i, j - 1)) / this->h), i, j);
		}
	}

	this->pressure->replaceDataVector(this->poissonSolver(100, 0.0000000001));

}

std::vector<std::vector<double>> Mesh::poissonSolver(int maxIterations, double epsilon)
{
	std::vector<std::vector<double>> residualVector;
	double r = 1 + epsilon;
	int counter = 0;

	std::vector<std::vector<double>> rightHandSide;
	std::vector<std::vector<double>> phi;

	phi = this->pressure->getData();
	rightHandSide = this->pressureRightHandSide->getData();

	while (counter < maxIterations && r > epsilon)
	{
		phi = multiGridSolver(phi, rightHandSide, this->h);
		residualVector = this->centerResidual(phi, rightHandSide, this->h);
		r = this->absoluteMax(residualVector);
		++counter;
	}
	return phi;
}

std::vector<std::vector<double>> Mesh::multiGridSolver(std::vector<std::vector<double>> &phi, std::vector<std::vector<double>> &rightHandSide, double h)
{
	double M = phi.size() - 2;
	double N = phi[0].size() - 2;
	std::vector<std::vector<double>> rightHandSideTemp;
	std::vector<std::vector<double>> phiTemp;

	phi = centerGS(phi, rightHandSide,h, 1);

	if (std::fmod(M, 2.0) == 0 && std::fmod(N, 2.0) == 0)
	{
		phiTemp.resize(M / 2 + 2);
		for (int i = 0; i < M / 2 + 2; ++i)
		{
			phiTemp[i].resize(N / 2 + 2);
		}

		rightHandSideTemp = this->centerResidual(phi, rightHandSide, h);
		rightHandSideTemp = this->centerRestrictCells(rightHandSideTemp);
		phiTemp = this->multiGridSolver(phiTemp, rightHandSideTemp, 2 * h);
		phiTemp = this->centerProlongCells(phiTemp);
		phi = this->addVectors(phi, phiTemp);
		this->setGSBoundaryConditions(phi);
		phi = this->centerGS(phi, rightHandSide, h, 1);
	}

	return phi;
}

std::vector<std::vector<double>> Mesh::centerGS(std::vector<std::vector<double>> &phi, std::vector<std::vector<double>> &rightHandSide, double h, int numIterations)
{
	double M = phi.size() - 2;
	double N = phi[0].size() - 2;
	this->setGSBoundaryConditions(phi);

	for (int k = 0; k < numIterations; ++k)
	{
		for (int j = 1; j < N + 1; ++j)
		{
			for (int i = 1; i < M + 1; ++i)
			{
				phi[i][j] = 0.25 * (phi[i - 1][j] + phi[i + 1][j] + phi[i][j - 1] + phi[i][j + 1]) - .25 * std::pow(h, 2) * rightHandSide[i][j];
			}
		}
		this->setGSBoundaryConditions(phi);
	}

	return phi;
}

std::vector<std::vector<double>> Mesh::centerResidual(std::vector<std::vector<double>> &phi, std::vector<std::vector<double>> &rightHandSide, double h)
{

	double M = phi.size() - 2;
	double N = phi[0].size() - 2;
	std::vector<std::vector<double>> residual;
	residual.resize(M + 2);
	for (int i = 0; i < M + 2; ++i)
	{
		residual[i].resize(N + 2);
	}

	for (int j = 1; j < N + 1; ++j)
	{
		for (int i = 1; i < M + 1; ++i)
		{
			residual[i][j] = rightHandSide[i][j] - ((phi[i - 1][j] - 2 * phi[i][j] + phi[i + 1][j]) / std::pow(h, 2) + (phi[i][j - 1] - 2 * phi[i][j] + phi[i][j + 1]) / std::pow(h, 2));
		}
	}

	return residual;
}

std::vector<std::vector<double>> Mesh::centerRestrictCells(std::vector<std::vector<double>> &vector)
{
	double M = vector.size() / 2 - 1;
	double N = vector[0].size() / 2 - 1;
	std::vector<std::vector<double>> restrictedVector;
	restrictedVector.resize(M + 2);
	for (int i = 0; i < M + 2; ++i)
	{
		restrictedVector[i].resize(N + 2);
	}

	for (int j = 1; j < N + 1; ++j)
	{
		for (int i = 1; i < M + 1; ++i)
		{
			restrictedVector[i][j] = 1.0 / 4.0 * (vector[2 * i - 2 + 1][2 * j - 2 + 1] + vector[2 * i - 1 + 1][2 * j - 2 + 1] + vector[2 * i - 2 + 1][2 * j - 1 + 1] + vector[2 * i - 1 + 1][2 * j - 1 + 1]);
		}
	}
	this->setGSBoundaryConditions(restrictedVector);

	return restrictedVector;
}

std::vector<std::vector<double>> Mesh::centerProlongCells(std::vector<std::vector<double>> &vector)
{
	double M = vector.size() - 2;
	double N = vector[0].size() - 2;
	std::vector<std::vector<double>> prolongedVector;
	prolongedVector.resize(2 * M + 2);
	for (int i = 0; i < 2 * M + 2; ++i)
	{
		prolongedVector[i].resize(2 * N + 2);
	}

	for (int j = 1; j < N + 1; ++j)
	{
		for (int i = 1; i < M + 1; ++i)
		{
			
			prolongedVector[2 * i - 2 + 1][2 * j - 2 + 1] = vector[i][j];
			prolongedVector[2 * i - 1 + 1][2 * j - 2 + 1] = vector[i][j];
			prolongedVector[2 * i - 2 + 1][2 * j - 1 + 1] = vector[i][j];
			prolongedVector[2 * i - 1 + 1][2 * j - 1 + 1] = vector[i][j];
		}
	}

	this->setGSBoundaryConditions(prolongedVector);

	return prolongedVector;
}

std::vector<std::vector<double>> Mesh::addVectors(std::vector<std::vector<double>> &vector1, std::vector<std::vector<double>> &vector2)
{
	if (vector1.size() != vector2.size() || vector1[0].size() != vector2[0].size()) 
	{
		throw std::invalid_argument("Vectors must be of the same size");
	}
	
	for (int j = 0; j < vector1.size(); ++j)
	{
		for (int i = 0; i < vector1[0].size(); ++i)
		{
			vector1[i][j] += vector2[i][j];
		}
	}
	return vector1;
}

void Mesh::setGSBoundaryConditions(std::vector<std::vector<double>> &phi)
{
	double M = phi.size() - 2;
	double N = phi[0].size() - 2;

	for (int j = 0; j < N + 2; ++j)
	{
		phi[M + 1][j] = phi[M][j];
		phi[0][j] = phi[1][j];
	}
	for (int i = 0; i < M + 2; ++i)
	{
		phi[i][0] = phi[i][1];
		phi[i][N + 1] = phi[i][N];
	}
}
