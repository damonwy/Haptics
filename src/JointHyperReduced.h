#pragma once

#ifndef REDUCEDCOORD_SRC_JOINTHYPERREDUCED_H_
#define REDUCEDCOORD_SRC_JOINTHYPERREDUCED_H_

#include "Joint.h"



class JointHyperReduced : public Joint {

public:
	JointHyperReduced() {}
	JointHyperReduced(std::shared_ptr<Body> body0, std::shared_ptr<Body> body1, std::shared_ptr<Joint> parent = nullptr);
	void load(const std::string &RESOURCE_DIR, std::string joint_shape);
	virtual ~JointHyperReduced() {}
	void update_();
protected:

	void draw_(std::shared_ptr<MatrixStack> MV,
		const std::shared_ptr<Program> prog,
		const std::shared_ptr<Program> prog2,
		std::shared_ptr<MatrixStack> P) const;
};

#endif // REDUCEDCOORD_SRC_JOINTHYPERREDUCED_H_