#pragma once

#ifndef REDUCEDCOORD_SRC_JOINTREVOLUTE_H_
#define REDUCEDCOORD_SRC_JOINTREVOLUTE_H_

#include "Joint.h"

class SE3;
class Body;

class JointRevolute : public Joint {

public:
	JointRevolute() {}
	JointRevolute(std::shared_ptr<Body> body, Eigen::Vector3d axis, std::shared_ptr<Joint> parent = nullptr);
	void load(const std::string &RESOURCE_DIR, std::string joint_shape);
	virtual ~JointRevolute() {}
	void update_();
protected:
	//virtual void computeHyperReducedJacobian_(Eigen::MatrixXd &JrR, Eigen::MatrixXd &JrR_select);

	void draw_(std::shared_ptr<MatrixStack> MV, 
		const std::shared_ptr<Program> prog, 
		const std::shared_ptr<Program> prog2, 
		std::shared_ptr<MatrixStack> P) const;
};

#endif // REDUCEDCOORD_SRC_JOINTREVOLUTE_H_