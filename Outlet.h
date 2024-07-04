#pragma once
#include "Opening.h"
class Outlet : public Opening
{
private:
	bool isOpen = true;
public:
	using Opening::Opening;
	bool getIsOpen() { return isOpen; }
	void toggleOpening();
};