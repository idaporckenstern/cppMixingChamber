#include <iostream>
#include <vector>
#include "Mesh.h"
#include "Functions.h"
#include "Opening.h"
#include "Outlet.h"
#include "WriteData.h"
#include <cmath>

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

    double plotTimes[6] = { 2, 5, 10, 15, 20, 24 };
    //std::vector<double> plotTimes = linspace(0.0, tEnd, tEnd * 10);
    //int size = plotTimes.size();
    //int size = 6;
    
    
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
    
    
    int plotCount = 1;
    
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

    //mesh.writeTimes(plotTimes);
    mesh.solveR();
    //mesh.writeData(0.0);


    while (mesh.getT() < mesh.getTEnd())
    {
        mesh.setDT(mesh.solveDT(CFL));
        if ((mesh.getT() < plotTimes[plotCount]) && (mesh.getT() + mesh.getDT() > plotTimes[plotCount]))
        {
            mesh.setDT(plotTimes[plotCount] - mesh.getT());
            ++plotCount;
            mesh.setDoPlot(true);
        }

        mesh.stepForward();

        //check if outlets are open or closed
        if (mesh.getT() < 10 || std::fmod(mesh.getT() - 10, 6) < 3)
        {
            mesh.setOpeningIsOpen(0, true);
        }
        else
        {
            mesh.setOpeningIsOpen(0, false);
        }

        if (mesh.getT() < 10 || std::fmod(mesh.getT() - 10, 6) >= 3)
        {
            mesh.setOpeningIsOpen(1, true);
        }
        else
        {
            mesh.setOpeningIsOpen(1, false);
        }
        //std::cout << "outlet 1: " << outlets[0].getIsOpen() << std::endl; 
        //std::cout << "outlet 2: " << outlets[1].getIsOpen() << std::endl;

        mesh.setBoundaryConditionsU(lambda, inletConditions, mesh.getT() + mesh.getDT());
        mesh.setBoundaryConditionsV(lambda, inletConditions, mesh.getT() + mesh.getDT());
        mesh.setBoundaryConditionsY(lambda, inletConditions, mesh.getT() + mesh.getDT());
        mesh.correctVelocities();
        mesh.setBoundaryConditionsU(lambda, inletConditions, mesh.getT());
        mesh.setBoundaryConditionsV(lambda, inletConditions, mesh.getT());

        if (mesh.getDoPlot())
        {
            mesh.solveR();
            //mesh.writeData(plotTimes[plotCount - 1]);
            mesh.setDoPlot(false);
        }

        mesh.setT(mesh.getT() + mesh.getDT());
        //mesh.testing();
        //std::cout << mesh.getT() << std::endl; 
        std::cout << mesh.getDT() << std::endl;

        
    }

    mesh.testing();
    
}


