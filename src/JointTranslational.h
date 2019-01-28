#pragma once

#ifndef REDUCEDCOORD_SRC_JOINTTRANSLATIONAL_H_
#define REDUCEDCOORD_SRC_JOINTTRANSLATIONAL_H_
#include "Joint.h"

#include "SE3.h"
#include "Shape.h"
#include "Body.h"

class JointTranslational : public Joint {

public:
	JointTranslational() {}
	JointTranslational(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr) :
		Joint(body, 3, parent)
	{

	}
	void load(const std::string &RESOURCE_DIR, std::string joint_shape) {
		//m_body->setJoint(getJoint());
		m_jointShape = std::make_shared<Shape>();
		m_jointShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	}

	virtual ~JointTranslational() {}

	void update_() {
		m_Q = Matrix4d::Identity();
		m_Q.block<3, 1>(0, 3) = m_q.topRows(3);
		m_S(3, 0) = 1.0;
		m_S(4, 1) = 1.0;
		m_S(5, 2) = 1.0;
	}

protected:

	void draw_(std::shared_ptr<MatrixStack> MV,
		const std::shared_ptr<Program> prog,
		const std::shared_ptr<Program> prog2,
		std::shared_ptr<MatrixStack> P) const {

	}



};



#endif // REDUCEDCOORD_SRC_JOINTTRANSLATIONAL_H_