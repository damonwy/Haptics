#pragma once
// SpringNull 

#ifndef REDUCEDCOORD_SRC_SPRINGNULL_H_
#define REDUCEDCOORD_SRC_SPRINGNULL_H_

#include "Spring.h"
class SpringNull : public Spring {

public:
	SpringNull() : Spring() {}
	virtual ~SpringNull() {}

};

#endif // REDUCEDCOORD_SRC_SPRINGNULL_H_