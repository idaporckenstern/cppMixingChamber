#include "WriteData.h"

void WriteData::Write(std::vector<std::vector<double>> data, double time, std::string dataType)
{
	this->time = std::to_string(time);
	this->dataFileName = "../Processed Data/Data/";
	this->dataFileName.append(dataType);
	this->dataFileName.append("t=");
	this->dataFileName.append(this->time);
	this->dataFileName.append("Data.txt");
	std::ofstream dataFile(this->dataFileName);

	for (int i = 0; i < data.size(); ++i)
	{
		for (int j = 0; j < data[0].size(); ++j)
		{
			dataFile << data[i][j];
			if (j < data[0].size() - 1)
			{
				dataFile << ",";
			}
		}
		dataFile << std::endl;
	}
}

void WriteData::WriteTimes(double times[], int size)
{
	std::ofstream timesFile("../Processed Data/Data/times.txt");

	for (int i = 0; i < size; ++i)
	{
		timesFile << std::setprecision(6) << std::fixed << times[i];
		if (i < size - 1)
		{
			timesFile << ",";
		}
	}

	timesFile.close();
}

void WriteData::WriteTimes(std::vector<double> times)
{
	std::ofstream timesFile("../Processed Data/Data/times.txt");

	for (int i = 0; i < times.size(); ++i)
	{
		timesFile << std::setprecision(6) << std::fixed << times[i];
		if (i < times.size() - 1)
		{
			timesFile << ",";
		}
	}

	timesFile.close();
}
