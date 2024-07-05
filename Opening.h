#pragma once
#include "Functions.h"

class Opening
{
protected:
	DataPoint startingPoint;
	DataPoint endPoint;

public:
	Opening(DataPoint startingPoint, DataPoint endPoint);
	DataPoint getStartingPoint();
	DataPoint getEndPoint();

};
