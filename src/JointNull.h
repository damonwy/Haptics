#pragma once

#ifndef REDUCEDCOORD_SRC_JOINTNULL_H_
#define REDUCEDCOORD_SRC_JOINTNULL_H_

#include "Joint.h"

class JointNull :public Joint {
public:
	JointNull() {}
	virtual ~JointNull() {}
	void update() {}
};

#endif // REDUCEDCOORD_SRC_JOINTNULL_H_