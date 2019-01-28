#pragma once

#ifndef REDUCEDCOORD_SRC_JOINTREVOLUTEHYPERREDUCED_H_
#define REDUCEDCOORD_SRC_JOINTREVOLUTEHTPERREDUCED_H_

#include "JointRevolute.h"

#include "Body.h"

class JointRevoluteHyperReduced : public JointRevolute {

public:
	JointRevoluteHyperReduced() {}
	JointRevoluteHyperReduced(std::shared_ptr<Body> body, Eigen::Vector3d axis, std::shared_ptr<Joint> friend_joint, double scalar, std::shared_ptr<Joint> parent = nullptr):
	JointRevolute(body, axis, parent), m_friend_joint(friend_joint), m_scalar(scalar)
	{
	}

	void countHRDofs(int &nR) {
	}

protected:
	void init_() {
		idxHR = m_friend_joint->idxHR;
	}

	std::shared_ptr<Joint> m_friend_joint;
	double m_scalar;

	void computeHyperReducedJacobian_(Eigen::MatrixXd &JrR, Eigen::MatrixXd &JrR_select) {
		JrR(idxR, idxHR) = m_scalar;
		JrR_select(idxR, idxHR) = 0.0;
	}
};

#endif // REDUCEDCOORD_SRC_JOINTREVOLUTEHYPERREDUCED_H_