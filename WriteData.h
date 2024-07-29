#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "Functions.h"

class WriteData
{
public:
	WriteData();


	void writeTimeJSON(double plotTimes[], int plotTimesSize, std::vector<double> animationTimes);
	void writeDataJSON(std::vector<std::vector<double>> data, FileType fileType);
	void closeFiles();

private:
	static inline int plotCounter = 0;
	static inline int animationCounter = 0;


	std::ofstream timeJSON;
	std::ofstream uPlotJSON;
	std::ofstream uAnimationJSON;
	std::ofstream vPlotJSON;
	std::ofstream vAnimationJSON;
	std::ofstream YPlotJSON;
	std::ofstream YAnimationJSON;
	std::ofstream RPlotJSON;
	std::ofstream RAnimationJSON;
};

