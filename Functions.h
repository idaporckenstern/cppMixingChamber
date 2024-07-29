#pragma once
#include <vector>
#include <string>
#include <iostream>

struct DataPoint
{
    double x;
    double y;

    DataPoint(double x, double y);
};

struct ServerAddress
{
    std::string ipAddress;
    int port;
};

enum class FileType
{
    uAnimation,
    uPlot,
    vAnimation,
    vPlot,
    YAnimation,
    YPlot,
    RAnimation,
    RPlot,
};

std::string fileTypeToString(FileType fileType);

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}
