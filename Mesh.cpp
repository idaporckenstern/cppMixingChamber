#include "Mesh.h"

Mesh::Mesh(double M, double N, double lengthX, double lengthY, double Re, double Sc, double t, double tEnd, std::vector<Opening> &inlets, std::vector<Outlet> &outlets)
	: M(M), N(N), lengthX(lengthX), lengthY(lengthY), h(lengthX/M), Re(Re), Sc(Sc), t(t), tEnd(tEnd),
	uVelocity(std::make_unique<DataVector>(linspace(0.0, lengthX, M + 1), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	vVelocity(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0, lengthY, N + 1))),
	Y(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	pressure(std::make_unique<DataVector>(linspace(0 - h / 2, lengthX + h / 2, M + 2), linspace(0.0 - h / 2.0, lengthY + h / 2.0, N + 2))),
	inlets(inlets), outlets(outlets)

{
	uVelocity->initDataVector();
	vVelocity->initDataVector();
	Y->initDataVector();
	pressure->initDataVector();
}

void Mesh::testing()
{

	std::cout << "U: \n";
	uVelocity->testing();
	std::cout << std::endl;
	std::cout << "V: \n";
	vVelocity->testing();
	std::cout << std::endl;
	std::cout << "Y: \n";
	Y->testing();
	std::cout << std::endl;
	std::cout << "pressure: \n";
	pressure->testing();
	std::cout << std::endl;
}