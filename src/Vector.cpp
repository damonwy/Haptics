#include "rmpch.h"
#include "Vector.h"
#include "Node.h"
#include "Body.h"

using namespace std;
using namespace Eigen;

Vector::Vector() :
	dir(0.0, 0.0, 0.0),
	dir0(0.0, 0.0, 0.0),
	fixed(false)
{
}

Vector::Vector(shared_ptr<Node> p, shared_ptr<Body> body, Vector3d dir):
m_p(p), m_body(body), dir0(dir), dir(dir)
{
}

void Vector::reset()
{
	dir = dir0;
}

void Vector::update(Matrix4d E) {
	Vector4d pos;
	pos.segment<3>(0) = this->dir0;
	pos(3) = 0.0;
	pos = E * pos;
	this->dir = pos.segment<3>(0);
}

void Vector::update() {
	if (m_body != nullptr) {
		Vector4d pos;
		pos.segment<3>(0) = this->dir0;
		pos(3) = 0.0;
		pos = m_body->E_wi * pos;
		this->dir = pos.segment<3>(0);
	}
}

void Vector::draw(shared_ptr<MatrixStack> MV, shared_ptr<MatrixStack> P, const shared_ptr<Program> prog) const
{
	if (m_p) {
		prog->bind();
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		MV->pushMatrix();
		glColor3f(0.8f, 0.7f, 0.0f);
		glLineWidth(3);
		glBegin(GL_LINES);
		Vector3f p0 = m_p->x.cast<float>();
		glVertex3f(p0(0), p0(1), p0(2));
		Vector3f p1 = (p0 + this->dir.cast<float>());
		glVertex3f(p1(0), p1(1), p1(2));
		glEnd();
		MV->popMatrix();
		prog->unbind();
	}
}
