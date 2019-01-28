#pragma once

#ifndef MUSCLEMASS_SRC_JOINTFIXED_H_
#define MUSCLEMASS_SRC_JOINTFIXED_H_
#define EIGEN_DONT_ALIGN_STATICALLY

#include <Eigen/Dense>
#include "Joint.h"
#include "MLCommon.h"

class Body;

class JointFixed : public Joint {

public:
	JointFixed() {}
	JointFixed(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr) :
		Joint(body, 0, parent)
	{
	}

	virtual ~JointFixed() {}

protected:


	void update_(){
		//E_pj = E_pj0;
		m_Q = Matrix4d::Identity();
	}
};



#endif MUSCLEMASS_SRC_JOINTFIXED_H_