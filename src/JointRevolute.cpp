#include "rmpch.h"
#include "JointRevolute.h"
#include "Body.h"
#include "ConstraintPrescJoint.h"

using namespace std;
using namespace Eigen;

JointRevolute::JointRevolute(std::shared_ptr<Body> body, Eigen::Vector3d axis, std::shared_ptr<Joint> parent):
Joint(body, 1, parent)
{
	m_axis = axis;
}

void JointRevolute::load(const std::string &RESOURCE_DIR, std::string joint_shape) {

	m_jointShape = make_shared<Shape>();
	m_jointShape->loadMesh(RESOURCE_DIR + "sphere2.obj");

}

void JointRevolute::update_() {
	Matrix3d R = SE3::aaToMat(m_axis, m_q(0));
	Matrix4d Q;
	Q.setIdentity();
	Q.block<3, 3>(0, 0) = R;
	m_Q = Q;
	m_S.block<3, 1>(0, 0) = m_axis;
}

void JointRevolute::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	prog->bind();

	
	if (m_jointShape) {
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		glUniform3f(prog->getUniform("lightPos1"), 66.0f, 25.0f, 25.0f);
		glUniform1f(prog->getUniform("intensity_1"), 0.6f);
		glUniform3f(prog->getUniform("lightPos2"), -66.0f, 25.0f, 25.0f);
		glUniform1f(prog->getUniform("intensity_2"), 0.2f);
		glUniform1f(prog->getUniform("s"), 300.0f);
		glUniform3f(prog->getUniform("ka"), 0.2f, 0.2f, 0.2f);
		glUniform3f(prog->getUniform("kd"), 0.8f, 0.7f, 0.7f);
		glUniform3f(prog->getUniform("ks"), 1.0f, 0.9f, 0.8f);

		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wj));
		double alpha = 1.0;
		if (presc != nullptr && presc->activeER) {
			alpha = 2.0;
		}

		MV->scale(static_cast<float>(alpha * m_draw_radius));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_jointShape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();
}
//
//void JointRevolute::computeHyperReducedJacobian_(MatrixXd &JrR, MatrixXd &JrR_select) {
//	JrR(idxR, idxHR) = 1.0;
//	JrR_select(idxR, idxHR) = 1.0;
//}

