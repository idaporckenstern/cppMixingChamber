#include <iostream>
#include <vector>
#include "Mesh.h"
#include "Functions.h"
#include "Opening.h"
#include "Outlet.h"
#include <cmath>
#include <cstdlib>

int main()
{
    double M = 16;
    double N = M;
    double lengthX = 4;
    double lengthY = 4;
    double Re = 40;
    double Sc = 2;
    double t = 0;
    double tEnd = 25;
    double CFL = 0.8;

    int plotTimes[6] = { 2, 5, 10, 15, 20, 24 };

    auto lambda = [](double z, double L)
        {
            return 6 * z / L * (1 - z / L);
        };
    std::vector<double> inletConditions({ 1.0 / std::sqrt(2), 1.0, 1.0 / std::sqrt(2), 1.0 });

    std::vector<Opening> inlets({ 
        Opening(DataPoint(0.0, lengthY / 16), DataPoint(0.0, 3 * lengthY / 16)), 
        Opening(DataPoint(13 * lengthX / 16, 0.0), DataPoint(15 * lengthX / 16, 0.0)),
        Opening(DataPoint(lengthX, 13 * lengthX / 16), DataPoint(lengthX, 15 * lengthX / 16)),
        Opening(DataPoint(lengthX / 16, lengthY), DataPoint(3 * lengthX / 16, lengthY))
    });
    std::vector<Outlet> outlets({
        Outlet(DataPoint(7 * lengthX / 16, 0.0), DataPoint(9 * lengthX / 16, 0.0)),
        Outlet(DataPoint(7 * lengthX / 16, lengthY), DataPoint(9 * lengthX / 16, lengthY))
    });

    Mesh mesh = Mesh(M, N, lengthX, lengthY, Re, Sc, t, tEnd, inlets, outlets);
    
    if (mesh.getT() < 10 || std::fmod(mesh.getT() - 10, 6) < 3)
    {
        outlets[0].setOpening(true);
    }
    else
    {
        outlets[0].setOpening(false);
    }

    if (mesh.getT() < 10 || std::fmod(mesh.getT() - 10, 6) >= 3)
    {
        outlets[1].setOpening(true);
    }
    else
    {
        outlets[1].setOpening(false);
    }

    mesh.setBoundaryConditionsU(lambda, inletConditions, mesh.getT());
    mesh.setBoundaryConditionsV(lambda, inletConditions, mesh.getT());
    mesh.setBoundaryConditionsY(lambda, inletConditions, mesh.getT());

    
    int plotCount = 0;
    while (mesh.getT() < mesh.getTEnd())
    {
        double dt = mesh.solveDT(CFL);
        if ((mesh.getT() < plotTimes[plotCount]) && (mesh.getT() + dt > plotTimes[plotCount]))
        {
            dt = plotTimes[plotCount] - mesh.getT();
            ++plotCount;
            mesh.setDoPlot(true);
        }

        mesh.stepForward();

        //check if outlets are open or closed
        if (mesh.getT() < 10 || std::fmod(mesh.getT() - 10, 6) < 3)
        {
            outlets[0].setOpening(true);
        }
        else
        {
            outlets[0].setOpening(false);
        }

        if (mesh.getT() < 10 || std::fmod(mesh.getT() - 10, 6) >= 3)
        {
            outlets[1].setOpening(true);
        }
        else
        {
            outlets[1].setOpening(false);
        }

        mesh.setBoundaryConditionsU(lambda, inletConditions, mesh.getT() + dt);
        mesh.setBoundaryConditionsV(lambda, inletConditions, mesh.getT() + dt);
        mesh.setBoundaryConditionsY(lambda, inletConditions, mesh.getT() + dt);
        mesh.correctVelocities();
        mesh.setBoundaryConditionsU(lambda, inletConditions, mesh.getT());
        mesh.setBoundaryConditionsV(lambda, inletConditions, mesh.getT());
        mesh.setT(mesh.getT() + dt);
        std::cout << mesh.getT() << std::endl;
        
    }

    mesh.testing();
    
    
}


