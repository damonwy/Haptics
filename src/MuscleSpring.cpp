#include "rmpch.h"
#include "MuscleSpring.h"
#include "Body.h"
#include "Node.h"
#include "World.h"
#include "Joint.h"
#include "SolverSparse.h"
#include "Solver.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;
#define EPSILON 1e-8
//#define DEBUG

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
	// Compute the sum of J_i.transpose() * f_i from each node
	for (int i = 0; i < m_n_nodes; ++i) {
		VectorXd Jfi = m_nodes[i]->m_J.transpose() * m_nodes[i]->m * grav; // m_n_bodies x 1
		f += Jfi;
	}
}

void MuscleSpring::computeForceDamping_(Vector3d grav, Eigen::VectorXd & f, Eigen::MatrixXd & D)
{
}

void MuscleSpring::computeForceDampingSparse_(Vector3d grav, Eigen::VectorXd & f, std::vector<T>& D_)
{
}

void MuscleSpring::computeEnergies_(Vector3d grav, Energy & ener)
{
	for (int i = 0; i < m_n_nodes; ++i) {
		ener.V = ener.V - m_nodes[i]->m * m_nodes[i]->x.transpose() * grav;
	}
}

void MuscleSpring::computeJMJ_(Eigen::MatrixXd &JMJ, shared_ptr<World> world)
{	
	// This function goes first
	VectorXd y0(2 * world->nr), yi;
	world->getJoint0()->gatherDofs(y0, world->nr);
	computeJ(y0, world);
	for (int i = 0; i < m_n_nodes; ++i) {
		JMJ += m_nodes[i]->m * m_nodes[i]->m_J.transpose() * m_nodes[i]->m_J;
	}
}

void MuscleSpring::computeJ(VectorXd q0, std::shared_ptr<World> world) {
	VectorXd yi;
	for (int i = 0; i < m_n_bodies; i++) {
		int idx = m_bodies[i]->getJoint()->idxR;
		forward(idx, q0, world);
		backward(idx, q0, world);

		for (int j = 0; j < m_n_nodes; ++j) {
			Vector3d diff = (m_nodes[j]->x_f - m_nodes[j]->x_b) / (2 * EPSILON);
			m_nodes[j]->m_J.col(idx) = diff; 
		}
	}
	recover(q0, world);
}

void MuscleSpring::recover(Eigen::VectorXd q0, std::shared_ptr<World> world) {
	// Restore the states
	world->getJoint0()->scatterDofs(q0, world->nr);
	world->getMuscle0()->update();
}


void MuscleSpring::forward(int i, Eigen::VectorXd q0, std::shared_ptr<World> world) {
	q0(i) += EPSILON;	
	world->getJoint0()->scatterDofs(q0, world->nr);
	world->getMuscle0()->update();
	for (int t = 0; t < m_n_nodes; ++t) {
		m_nodes[t]->saveForwardPosition();
	}
}

void MuscleSpring::backward(int i, Eigen::VectorXd q0, std::shared_ptr<World> world) {
	q0(i) -= EPSILON;
	world->getJoint0()->scatterDofs(q0, world->nr);
	world->getMuscle0()->update();
	for (int t = 0; t < m_n_nodes; ++t) {
		m_nodes[t]->saveBackwardPosition();
	}
}

void MuscleSpring::computeJMJSparse_(std::vector<T>& J_)
{
}

void MuscleSpring::computeJMJdotqdot_(Eigen::VectorXd & f, const Eigen::VectorXd & qdot, shared_ptr<World> world, std::shared_ptr<SolverSparse> solver)
{
	// Get current states
	VectorXd q0 = solver->q0;
	VectorXd qdot0 = solver->qdot0;
	VectorXd y0(2*world->nr), yi, yj;
	y0 << q0, qdot0;

	// Compute Jdot = dJdq * qdot
	for (int i = 0; i < m_n_bodies; ++i) {
		// Compute Jdot for each node vector<3x2> 2 	
		int idx = m_bodies[i]->getJoint()->idxR;
		computedJdq(y0, world, true, idx);
		computedJdq(y0, world, false, idx);
		// Now m_dJdq 3 x dof x dof stores the Jacobian after perturbing each body a little bit
	}

	recover(y0, world);

	//Compute dJdq and Jdot for each node
	for (int t = 0; t < m_n_nodes; ++t) {
		auto node = m_nodes[t];
		node->m_Jdot.setZero();

		for (int s = 0; s < m_n_bodies; ++s) { // each slice of dJdq
			if (t == 1) {
				//cout << "node " << t << "slice " << s << endl << m_nodes[t]->m_dJdq[s] << endl << endl;
				//cout << "node " << t << "slice " << s << endl << m_nodes[t]->m_dJdq_b[s] << endl << endl;
			}

			node->m_dJdq[s] = (node->m_dJdq[s] - node->m_dJdq_b[s])/(2 *EPSILON);

#ifdef DEBUG			
			if (t == 1) {
				cout << "node " << t << "slice " << s << endl << m_nodes[t]->m_dJdq[s] << endl << endl;
		
			}
#endif
			node->m_Jdot += node->m_dJdq[s] * qdot0(s);
		}
	}

	for (int j = 0; j < m_n_nodes; ++j) {
		auto node = m_nodes[j];
		VectorXd fvec = node->m * node->m_J.transpose() * Matrix3d::Identity() * node->m_Jdot * qdot0; // need to check the sign later
		f -= fvec;	
#ifdef DEBUG
		if (j == 1) {
			cout << "f_" << j << endl << fvec << endl << endl;
			cout << "J " << endl << node->m_J << endl << endl;
			cout << "Jdot " << endl << node->m_Jdot << endl << endl;
			cout << "qdot " << endl << qdot0 << endl << endl;
		}
#endif
	}
	
#ifdef DEBUG
	Matrix2d testJ, testJdot;
	vector<Matrix2d> testdJdq;
	test(y0, testJ, testdJdq, testJdot);
	cout << "test: " << endl;
	cout << "J: " << endl << testJ << endl;
	cout << "Jdot: " << endl << testJdot << endl;
	cout << "dJdq 0: " << endl << testdJdq[0] << endl;
	cout << "dJdq 1: " << endl << testdJdq[1] << endl;
#endif

}

void MuscleSpring::test(Eigen::VectorXd q0, Eigen::Matrix2d &J, std::vector<Eigen::Matrix2d> &dJdq, Eigen::Matrix2d &Jdot) {
	// Use this to test the J, Jdot, dJdq of a material point in a spring 
	// Compare this with the fd results
	// Use Epsilon = 1e-5

	double s = 0.5; // point position in a spring
	double l = 10.0;
	double r = 5.0;
	q0(1) += M_PI / 2.0;

	// 
	double s0 = sin(q0(1));
	double c0 = cos(q0(1));
	double s01 = sin(q0(0) + q0(1));
	double c01 = cos(q0(0) + q0(1));

	// Unpack the data
	Vector2d q, qdot;
	qdot(0) = q0(3); 
	qdot(1) = q0(2);
	q(0) = q0(1); 
	q(1) = q0(0);

	J.resize(2, 2);
	J << -l * s0 - r * s01, -r * s01, l*c0 + r * c01, r*c01;
	J *= s;

	Matrix2d M;
	M = J.transpose() * J;

	dJdq.clear();
	Matrix2d Jdot_;
	Jdot_ << -l * c0 - r * c01, -r * c01, -l * s0 - r * s01, -r * s01;
	Jdot_ *= s;
	dJdq.push_back(Jdot_);

	Jdot_ << -r * c01, -r * c01, -r * s01, -r * s01; 
	Jdot_ *= s;
	dJdq.push_back(Jdot_);

	Jdot.resize(2, 2);
	Jdot.setZero();

	for (int i = 0; i < 2; ++i) {
		Jdot += dJdq[i] * qdot(i);
	}

}

void MuscleSpring::computedJdq(Eigen::VectorXd y0, std::shared_ptr<World> world, bool isForward, int idx) {
	if (isForward) {
		y0(idx) += EPSILON;
	}
	else {
		y0(idx) -= EPSILON;
	}

	for (int j = 0; j < m_n_bodies; j++) {
		int index = m_bodies[j]->getJoint()->idxR;
		forward(index, y0, world);
		backward(index, y0, world);

		for (int k = 0; k < m_n_nodes; ++k) {
			Vector3d diff = (m_nodes[k]->x_f - m_nodes[k]->x_b) / (2 * EPSILON);
			//cout << "layer: idx " << idx << "col" << index << endl << diff << endl;
			if (isForward) {
				m_nodes[k]->m_dJdq[idx].col(index) = diff; // attention: the corresponding rigid body
			}
			else {
				m_nodes[k]->m_dJdq_b[idx].col(index) = diff; // attention: the corresponding rigid body
			}
		}
	}
}
