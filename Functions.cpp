#include "Functions.h"

DataPoint::DataPoint(double x, double y) : x(x), y(y) {}

std::string fileTypeToString(FileType fileType) throw()
{
	switch (fileType)
	{
	case FileType::uPlot: return "\"u Plot\"";
	case FileType::uAnimation: return "\"u Animation\"";
	case FileType::vPlot: return "\"v Plot";
	case FileType::vAnimation: return "\"v Animation\"";
	case FileType::YPlot: return "\"Y Plot";
	case FileType::YAnimation: return "\"Y Animation\"";
	case FileType::RPlot: return "\"R Plot";
	case FileType::RAnimation: return "\"R Animation\"";
	default: throw std::invalid_argument("\"Unimplemented item\"");
	}
}
