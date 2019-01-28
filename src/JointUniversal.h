#pragma once
#include "Joint.h"
#ifndef REDUCEDCOORD_SRC_JOINTUNIVERSAL_H_
#define REDUCEDCOORD_SRC_JOINTUNIVERSAL_H_

#include "Body.h"
#include "SE3.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"
#include "ConstraintPrescJoint.h"


class JointUniversal : public Joint {

public:
	JointUniversal() {}
	JointUniversal(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr):
	Joint(body, 2, parent)
	{


	}

	virtual void load(const std::string &RESOURCE_DIR, std::string joint_shape) {
		m_jointShape = std::make_shared<Shape>();
		m_jointShape->loadMesh(RESOURCE_DIR + "sphere2.obj");

	}

	virtual ~JointUniversal() {}

protected:

	virtual void draw_(std::shared_ptr<MatrixStack> MV,
		const std::shared_ptr<Program> prog,
		const std::shared_ptr<Program> prog2,
		std::shared_ptr<MatrixStack> P) const {

		prog->bind();

		//if (m_jointShape) {
		//	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		//	glUniform3f(prog->getUniform("lightPos1"), 66.0f, 25.0f, 25.0f);
		//	glUniform1f(prog->getUniform("intensity_1"), 0.6f);
		//	glUniform3f(prog->getUniform("lightPos2"), -66.0f, 25.0f, 25.0f);
		//	glUniform1f(prog->getUniform("intensity_2"), 0.2f);
		//	glUniform1f(prog->getUniform("s"), 300.0f);
		//	glUniform3f(prog->getUniform("ka"), 0.2f, 0.2f, 0.2f);
		//	glUniform3f(prog->getUniform("kd"), 0.8f, 0.7f, 0.7f);
		//	glUniform3f(prog->getUniform("ks"), 1.0f, 0.9f, 0.8f);

		//	MV->pushMatrix();
		//	MV->multMatrix(eigen_to_glm(E_wj));
		//	//std::cout << E_wj << std::endl;
		//	double alpha = 1.0;

		//	if (presc!= nullptr && presc->activeER) {
		//		alpha = 2.0;
		//	}

		//	MV->scale(alpha * m_draw_radius);
		//	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		//	m_jointShape->draw(prog);
		//	MV->popMatrix();
		//}

		prog->unbind();
	}
};

#endif // REDUCEDCOORD_SRC_JOINTUNIVERAL_H_
