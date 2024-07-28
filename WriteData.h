#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

class WriteData
{
public:
	void Write(std::vector<std::vector<double>> data, double time, std::string dataType);
	void WriteTimes(double times[], int size);
	void WriteTimes(std::vector<double> times);

private:
	std::string time;
	std::string dataFileName;
};

