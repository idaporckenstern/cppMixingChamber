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
	std::cout << "Y: \n";
	Y->testing();
	std::cout << std::endl;
	//std::cout << "Y2: \n";
	//Y2->testing();
	std::cout << std::endl;
	std::cout << "pressure: \n";
	pressure->testing();
	std::cout << std::endl;

	
}

double Mesh::getDT(double CFL)
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
			std::cout << "\n";
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

double Mesh::sumX(std::vector<std::vector<double>>& data, int column, int start, int end)
{
	double sum = 0;
	for (int i = start; i < end; ++i)
	{
		sum += data[i][column];
	}
	return sum;
}

double Mesh::sumY(std::vector<std::vector<double>>& data, int row, int start, int end)
{
	double sum = 0;
	for (int j = start; j < end; ++j)
	{
		sum += data[row][j];
	}

	return sum;
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
}

void Mesh::conservationMassCorrection()
{

}
