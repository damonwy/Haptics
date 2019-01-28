#include "rmpch.h"
#include "CompSphere.h"

#include "Body.h"
#include "Node.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

CompSphere::CompSphere() {
	m_O = make_shared<Node>();
	m_O->x0.setZero();

}

CompSphere::CompSphere(std::shared_ptr<Body> parent, double r) :m_parent(parent), m_r(r){
	m_O = make_shared<Node>();
	m_O->x0.setZero();

}

void CompSphere::init() {
	m_shape->init();
	m_O->init();
}

void CompSphere::load(const string &RESOURCE_DIR) {
	m_shape = make_shared<Shape>();
	m_shape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	m_O->load(RESOURCE_DIR);
}


void CompSphere::update_() {
	E_wi = m_parent->E_wi * E_ji;
	m_O->update(E_wi);

}

void CompSphere::setTransform(Eigen::Matrix4d E) {

	E_ji = E;

}

void CompSphere::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P)const {

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
		m_O->draw(MV, prog);
		
		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wi));
		MV->scale(float(m_r));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();

}

