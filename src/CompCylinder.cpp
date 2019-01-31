#include "rmpch.h"
#include "CompCylinder.h"

#include "Body.h"
#include "Vector.h"
#include "Node.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

CompCylinder::CompCylinder() {
}

CompCylinder::CompCylinder(shared_ptr<Body> parent, double r) :  m_r(r), m_parent(parent){
}

void CompCylinder::init() {
	m_shape->init();
}

void CompCylinder::update_() {
	E_wi = m_parent->E_wi * E_ji;

	m_Z->update(E_wi);
	m_O->update(E_wi);

}

void CompCylinder::setTransform(Matrix4d E) {
	E_ji = E;
	E_wi = m_parent->E_wi * E_ji;
	
}

void CompCylinder::load(const string &RESOURCE_DIR, string shape) {
	m_shape = make_shared<Shape>();
	m_shape->loadMesh(RESOURCE_DIR + shape);
}

void CompCylinder::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P)const {
	prog->bind();
	if (m_shape) {
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
		MV->multMatrix(eigen_to_glm(E_wi));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shape->draw(prog);
		MV->popMatrix();

	}
	prog->unbind();

}
