#pragma once

#ifndef REDUCEDCOORD_SRC_JOINTFREE_H_
#define REDUCEDCOORD_SRC_JOINTFREE_H_

#include "Joint.h"

#include "SE3.h"
#include "Shape.h"
#include "Body.h"
#include "JointTranslational.h"
#include "JointSphericalExp.h"
#include "MatrixStack.h"
#include "Program.h"

class JointFree : public Joint {

public:
	JointFree() {}
	JointFree(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr):
	Joint(body, 6, parent)
	{
		m_jointS = std::make_shared<JointSphericalExp>(body);
		m_jointT = std::make_shared<JointTranslational>(body);

	}
	void load(const std::string &RESOURCE_DIR, std::string joint_shape) {
		m_body->setJoint(getJoint());
		m_jointShape = std::make_shared<Shape>();
		m_jointShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	}

	virtual ~JointFree() {}
	inline void reparam_() {
		m_jointS->reparam_();
		m_q.topRows(3) = m_jointS->m_q;
		m_qdot.topRows(3) = m_jointS->m_qdot;
	}	
	
	void update_() {
		m_jointS->m_q = m_q.topRows(3);	
		m_jointT->m_q = m_q.bottomRows(3);		
		m_jointS->m_qdot = m_qdot.topRows(3);
		m_jointT->m_qdot = m_qdot.bottomRows(3);
		//ok

		m_jointS->update_();
		m_jointT->update_();

		Matrix4d Q1 = m_jointS->m_Q;
		Matrix4d Q2 = m_jointT->m_Q;
		m_Q = Q1 * Q2;
		
		Eigen::MatrixXd S1 = m_jointS->m_S;
		Eigen::MatrixXd S2 = m_jointT->m_S;
		Matrix3d S1w = S1.block<3, 3>(0, 0);
		Vector3d p = Q2.block<3, 1>(0, 3);
		Matrix3d pbrac = SE3::bracket3(p);
		m_S.block<3, 3>(0, 0) = S1w;
		m_S.block<3, 3>(3, 0) = -pbrac * S1w;
		m_S.block<6, 3>(0, 3) = S2;
		///std::cout << "m_S" << std::endl << m_S << std::endl << std::endl;

		Eigen::MatrixXd dS1 = m_jointS->m_Sdot;
		Eigen::MatrixXd dS2 = m_jointT->m_Sdot;

		Matrix3d dS1w = dS1.block<3, 3>(0, 0);
		m_Sdot.block<3, 3>(0, 0) = dS1w;
		m_Sdot.block<3, 3>(3, 0) = -(SE3::bracket3(m_qdot.bottomRows(3))*S1w + pbrac * dS1w);
		m_Sdot.block<6, 3>(0, 3) = dS2;
		///std::cout << "m_Sdot" << std::endl << m_Sdot << std::endl << std::endl;

	}

protected:
	std::shared_ptr<JointSphericalExp> m_jointS;
	std::shared_ptr<JointTranslational> m_jointT;

	void draw_(std::shared_ptr<MatrixStack> MV,
		const std::shared_ptr<Program> prog,
		const std::shared_ptr<Program> prog2,
		std::shared_ptr<MatrixStack> P) const {
	}
};



#endif // REDUCEDCOORD_SRC_JOINTFREE_H_