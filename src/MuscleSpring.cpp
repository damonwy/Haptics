#include "rmpch.h"
#include "MuscleSpring.h"
#include "Body.h"
#include "Node.h"
#include "World.h"
#include "Joint.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;
#define EPSILON 1e-8

MuscleSpring::MuscleSpring(std::vector<std::shared_ptr<Body>> bodies, int n_nodes):
	Muscle(bodies, n_nodes)
{
}

void MuscleSpring::setMass(double mass)
{
	m_mass = mass; 
	double mass_i = m_mass / m_n_nodes;
	for (int i = 0; i < m_n_nodes; ++i) {
		m_nodes[i]->m = mass_i;
	}
}

void MuscleSpring::setAttachments(std::shared_ptr<Body> body0, Vector3d r0, std::shared_ptr<Body> body1, Vector3d r1)
{
	m_body0 = body0;
	m_body1 = body1;
	m_r0 = r0;
	m_r1 = r1;
}

void MuscleSpring::init_()
{
	for (int i = 0; i < m_n_nodes; i++) {
		m_nodes[i]->init();
	}

	update_();
}

void MuscleSpring::update_() {
	Matrix4d E0, E1;
	if (m_body0 == nullptr) {
		E0 = Matrix4d::Identity();
	}
	else {
		E0 = m_body0->E_wi;
	}

	if (m_body1 == nullptr) {
		E1 = Matrix4d::Identity();
	}
	else {
		E1 = m_body1->E_wi;
	}
	Vector4d r0, r1;
	r0.segment<3>(0) = m_r0;
	r1.segment<3>(0) = m_r1;
	r0(3) = 1.0;
	r1(3) = 1.0;
	Vector4d x0 = E0 * r0;
	Vector4d x1 = E1 * r1;

	// Set the nodal positions

	for (int i = 0; i < m_n_nodes; i++) {
		double s = double(i) / (m_n_nodes - 1);
		m_nodes[i]->x = (1 - s) * x0.segment<3>(0) + s * x1.segment<3>(0);
	}
}

void MuscleSpring::draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const
{
	// Draw nodes
	prog->bind();

	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	for (int i = 0; i < m_n_nodes; i++) {
		m_nodes[i]->draw(MV, prog);
	}

	MV->popMatrix();
	prog->unbind();

	// Draw line segments
	progSimple->bind();
	glUniformMatrix4fv(progSimple->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	glLineWidth(5);
	glBegin(GL_LINES);
	glColor3f(0.0f, 0.0f, 0.0f);

	for (int i = 0; i < m_n_nodes - 1; i++) {

		Vector3f x0 = m_nodes[i]->x.cast<float>();
		Vector3f x1 = m_nodes[i + 1]->x.cast<float>();

		glVertex3f(x0(0), x0(1), x0(2));
		glVertex3f(x1(0), x1(1), x1(2));
	}
	glEnd();
	progSimple->unbind();
}

void MuscleSpring::computeMass_(Eigen::MatrixXd & M)
{
}

void MuscleSpring::computeMassSparse_(std::vector<T>& M_)
{
}

void MuscleSpring::computeForce_(Vector3d grav, Eigen::VectorXd & f)
{
}

void MuscleSpring::computeForceDamping_(Vector3d grav, Eigen::VectorXd & f, Eigen::MatrixXd & D)
{
}

void MuscleSpring::computeForceDampingSparse_(Vector3d grav, Eigen::VectorXd & f, std::vector<T>& D_)
{
}

void MuscleSpring::computeEnergies_(Vector3d grav, Energy & ener)
{
}

void MuscleSpring::computeJMJ_(Eigen::MatrixXd &JMJ, shared_ptr<World> world)
{	
	// This function goes first

	VectorXd y0(2 * world->nr), yi;
	world->getJoint0()->gatherDofs(y0, world->nr);
	for (int i = 0; i < m_n_nodes; ++i) {
		m_nodes[i]->savePosition();
	}

	
	for (int i = 0; i < m_n_bodies; i++) {
		// For each related body, perturb 
		yi = y0;
		int idx = m_bodies[i]->getJoint()->idxR;
		yi(idx) += EPSILON;
		world->getJoint0()->scatterDofs(yi, world->nr);
		world->getMuscle0()->update();

		for (int j = 0; j < m_n_nodes; ++j) {
			Vector3d diff = (m_nodes[j]->x - m_nodes[j]->x_s) / EPSILON;
			m_nodes[j]->m_J.col(i) = diff;
		}
	}

	// Ji computed 
	// Restore the states
	world->getJoint0()->scatterDofs(y0, world->nr);
	world->getMuscle0()->update();

	Matrix3d I3 = Matrix3d::Identity();
	for (int i = 0; i < m_n_nodes; ++i) {
		JMJ += m_nodes[i]->m * m_nodes[i]->m_J.transpose() * I3 * m_nodes[i]->m_J;
	}

}

void MuscleSpring::computeJMJSparse_(std::vector<T>& J_)
{
}
