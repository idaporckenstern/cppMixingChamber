#include "Opening.h"
#include <iostream>

Opening::Opening(DataPoint startingPoint, DataPoint endPoint)
{
	this->startingPoint = startingPoint;
	this->endPoint = endPoint;
}

DataPoint Opening::getStartingPoint()
{
	return this->startingPoint;
}


DataPoint Opening::getEndPoint()
{
	return this->endPoint;
}
