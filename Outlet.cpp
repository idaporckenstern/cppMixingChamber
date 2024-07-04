#include "Outlet.h"


void Outlet::toggleOpening()
{
	if (this->isOpen)
	{
		isOpen = false;
	}
	else
	{
		isOpen = true;
	}
}


