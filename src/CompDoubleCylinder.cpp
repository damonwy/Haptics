#include "rmpch.h"
#include "CompDoubleCylinder.h"
#include "Body.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

CompDoubleCylinder::CompDoubleCylinder(): Comp() {
}

CompDoubleCylinder::CompDoubleCylinder(shared_ptr<Body> parentA, double rA, shared_ptr<Body> parentB, double rB) : 
Comp(), m_rA(rA), m_rB(rB), m_parentA(parentA), m_parentB(parentB)
{
	m_OA = make_shared<Node>();
	m_OB = make_shared<Node>();
}

void CompDoubleCylinder::load(const string &RESOURCE_DIR, string shapeA, string shapeB) {
	m_shapeA = make_shared<Shape>();
	m_shapeA->loadMesh(RESOURCE_DIR + shapeA);

	m_shapeB = make_shared<Shape>();
	m_shapeB->loadMesh(RESOURCE_DIR + shapeB);

	m_OA->r = 0.1;
	m_OB->r = 0.1;
	m_OA->load(RESOURCE_DIR);
	m_OB->load(RESOURCE_DIR);
}

void CompDoubleCylinder::init() {
	m_shapeA->init();
	m_shapeB->init();
	m_OA->init();
	m_OB->init();
}

void CompDoubleCylinder::update_() {
	E_wiA = m_parentA->E_wi * E_jiA;
	E_wiB = m_parentB->E_wi * E_jiB;
	
	m_OA->update(E_wiA);
	m_OB->update(E_wiB);
	m_ZA->update(E_wiA);
	m_ZB->update(E_wiB);
}

void CompDoubleCylinder::setTransformA(Matrix4d E) {
	E_jiA = E;
}

void CompDoubleCylinder::setTransformB(Matrix4d E) {

	E_jiB = E;
}

void CompDoubleCylinder::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P)const {
	prog->bind();
	if (m_shapeA && m_shapeB) {
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
		MV->multMatrix(eigen_to_glm(E_wiA));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shapeA->draw(prog);
		MV->popMatrix();

		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wiB));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_shapeB->draw(prog);
		MV->popMatrix();

		m_OA->draw(MV, prog);
		m_OB->draw(MV, prog);

	}
	prog->unbind();
}
