#pragma once

#include "Constraint.h"

class ConstraintNull : public Constraint
{
public:
	ConstraintNull() : Constraint(0, 0, 0, 0) {

		m_name = "null";
	}

	virtual ~ConstraintNull() {}
};