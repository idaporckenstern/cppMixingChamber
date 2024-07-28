#include "Opening.h"
#include <iostream>

Opening::Opening(DataPoint startingPoint, DataPoint endPoint) : startingPoint(startingPoint), endPoint(endPoint)
{
	
}

DataPoint Opening::getStartingPoint()
{
	return this->startingPoint;
}


DataPoint Opening::getEndPoint()
{
	return this->endPoint;
}
