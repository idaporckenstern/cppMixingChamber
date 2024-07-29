#include "WriteData.h"

WriteData::WriteData(): timeJSON("../Processed Data/Data/times.json"),uPlotJSON("../Processed Data/Data/u_plot_data.json"),uAnimationJSON("../Processed Data/Data/u_animation_data.json"),
	vPlotJSON("../Processed Data/Data/v_plot_data.json"), vAnimationJSON("../Processed Data/Data/v_animation_data.json"),
	YPlotJSON("../Processed Data/Data/u_plot_data.json"), YAnimationJSON("../Processed Data/Data/Y_animation_data.json"),
	RPlotJSON("../Processed Data/Data/u_plot_data.json"), RAnimationJSON("../Processed Data/Data/R_animation_data.json")
{
	timeJSON << "{\n";
	uPlotJSON << "{\n";
	uAnimationJSON << "{\n";
}


void WriteData::writeTimeJSON(double plotTimes[], int plotTimesSize, std::vector<double> animationTimes)
{
	this->timeJSON << "\"plot times\": [";

	for (int i = 0; i < plotTimesSize; ++i)
	{
		this->timeJSON << std::setprecision(6) << std::fixed << plotTimes[i];
		if (i != plotTimesSize - 1)
		{
			this->timeJSON << ",";
		}
	}
	this->timeJSON << "],\n \"animation times\": [";
	
	for (int i = 0; i < animationTimes.size(); ++i)
	{
		this->timeJSON << std::setprecision(6) << std::fixed << animationTimes[i];
		if (i != animationTimes.size() - 1)
		{
			this->timeJSON << ",";
		}
	}
	this->timeJSON << "]";
}

void WriteData::writeDataJSON(std::vector<std::vector<double>> data, FileType fileType)
{
	std::ofstream* file = nullptr;

	switch (fileType)
	{
	case FileType::uPlot: file = &this->uPlotJSON;
		break;
	case FileType::uAnimation: file = &this->uAnimationJSON;
		break;
	case FileType::vPlot: file = &this->vPlotJSON;
		break;
	case FileType::vAnimation: file = &this->vAnimationJSON;
		break;
	case FileType::YPlot: file = &this->YPlotJSON;
		break;
	case FileType::YAnimation: file = &this->YAnimationJSON;
		break;
	case FileType::RPlot: file = &this->RPlotJSON;
		break;
	case FileType::RAnimation: file = &this->RAnimationJSON;
		break;
	}

	if (
		fileType == FileType::uPlot ||
		fileType == FileType::vPlot ||
		fileType == FileType::YPlot ||
		fileType == FileType::RPlot
		)
	{
		*file << "\"" << this->plotCounter << "\":[\n";
	}
	else
	{
		*file << "\"" << this->animationCounter << "\":[\n";
	}
	
	if (fileType == FileType::RPlot)
	{
		plotCounter++;
	}
	if (fileType == FileType::RAnimation)
	{
		animationCounter++;
	}

	for (int i = 0; i < data.size(); ++i)
	{
		*file << "[";
		for (int j = 0; j < data[0].size(); ++j)
		{
			*file << data[i][j];
			if (j < data[0].size() - 1)
			{
				*file << ",";
			}
		}
		*file << "]";
		if (i < data.size() - 1)
		{
			*file << ",\n";
		}
		else
		{
			*file << "\n";
		}
	}
	*file << "],\n";
	
}


void WriteData::closeFiles()
{
	this->timeJSON << "\n}";

	long pos = this->uPlotJSON.tellp();
	this->uPlotJSON.seekp(pos - 1);
	this->uPlotJSON.write("", 1);
	this->uPlotJSON << "\n}";

	 pos = this->uAnimationJSON.tellp();
	this->uAnimationJSON.seekp(pos - 3);
	this->uAnimationJSON.write(" ", 1);
	this->uAnimationJSON << "\n}";

	this->timeJSON.close();
	this->uPlotJSON.close();
	this->uAnimationJSON.close();
	this->vPlotJSON.close();
	this->vAnimationJSON.close();
	this->YPlotJSON.close();
	this->YAnimationJSON.close();
	this->RPlotJSON.close();
	this->RAnimationJSON.close();
}






