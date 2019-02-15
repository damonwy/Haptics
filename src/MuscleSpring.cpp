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
		//cout << "yi" << yi << endl;
		world->getJoint0()->scatterDofs(yi, world->nr);
		world->getMuscle0()->update();

		for (int j = 0; j < m_n_nodes; ++j) {
			Vector3d diff = (m_nodes[j]->x - m_nodes[j]->x_s) / EPSILON;
			m_nodes[j]->m_J.col(idx) = diff; // attention: the corresponding rigid body
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

void MuscleSpring::computeJMJdotqdot_(Eigen::VectorXd & f, const Eigen::VectorXd & qdot, shared_ptr<World> world)
{
	for (int i = 0; i < m_n_nodes; ++i) {	
		auto node = m_nodes[i];
		node->savePosition();
		//node->m_Jdot.clear();
	}

	VectorXd y0(2 * world->nr), ydot0(2 * world->nr), yi, yj;
	world->getJoint0()->gatherDofs(y0, world->nr);
	world->getJoint0()->gatherDDofs(ydot0, world->nr);

	for (int i = 0; i < m_n_bodies; ++i) {
		// Compute Jdot for each node vector<3x2> 2 
		yi = y0;
		int idx = m_bodies[i]->getJoint()->idxR;
		yi(idx) += EPSILON;
		// before perturb
		world->getJoint0()->scatterDofs(yi, world->nr);
		world->getMuscle0()->update();
		for (int t = 0; t < m_n_nodes; ++t) {
			m_nodes[t]->savePositionForJdot();
		}
		// 

		for (int j = 0; j < m_n_bodies; j++) {
			// For each related body, perturb 
			yj = yi;
			int index = m_bodies[j]->getJoint()->idxR;
			yj(index) += EPSILON;
			//cout << "yj" << yj << endl;
			world->getJoint0()->scatterDofs(yj, world->nr);
			world->getMuscle0()->update();

			for (int k = 0; k < m_n_nodes; ++k) {
				Vector3d diff = (m_nodes[k]->x - m_nodes[k]->x_ss) / EPSILON;
				m_nodes[k]->m_Jdot[idx].col(index) = diff;
				
				 // attention: the corresponding rigid body
			}
		}
	}

	// Every node has m_dJdq now
	VectorXd f_qvv0 = f;
	VectorXd f_qvv1 = f;

	for (int i = 0; i < m_n_bodies; ++i) {
		int idx = m_bodies[i]->getJoint()->idxR;
		for (int j = 0; j < m_n_nodes; ++j) {
			auto node = m_nodes[j];
			MatrixXd Mnew = node->m * node->m_Jdot[idx].transpose() * node->m_Jdot[idx];
			MatrixXd Mold = node->m * node->m_J.transpose() * node->m_J;
			MatrixXd diff_M = (Mnew - Mold) / EPSILON;
			node->m_dMdq[idx] = diff_M;

			if (idx == 1) {
				node->m_dMdq[idx].setZero();
			}			
			cout << "dMdq_" << idx << endl << node->m_dMdq[idx] << endl << endl;

			f_qvv0(idx) += 0.5* qdot.transpose() * node->m_dMdq[idx] * qdot;
		}
	}
	cout << "fqvv0:" << endl << f_qvv0 << endl << endl;

	Matrix2d dIdtheta;
	dIdtheta.setZero();

	for (int i = 0; i < m_n_bodies; ++i) {
		int idx = m_bodies[i]->getJoint()->idxR;
		for (int j = 0; j < m_n_nodes; ++j) {
			auto node = m_nodes[j];
			dIdtheta += node->m_dMdq[idx] * qdot(idx);
			//cout << "idx: " << idx << "j: " << j << "dIdtheta " << endl <<node->m_dMdq[idx] * qdot(idx) << endl;
			f_qvv1 += node->m_dMdq[idx]*qdot(idx)*qdot;
		}
	}
	cout << "dIdtheta: " << endl << dIdtheta << endl;
	cout << "fqvv1:" << endl << f_qvv1 << endl << endl;

	//VectorXd y0(2 * world->nr), ydot0(2*world->nr), yi, yj;
	//world->getJoint0()->gatherDofs(y0, world->nr);
	//world->getJoint0()->gatherDDofs(ydot0, world->nr);

	////cout << "ydot0" << endl << ydot0 << endl;
	////cout << "qdot0" << qdot << endl;


	//for (int i = 0; i < m_n_nodes; ++i) {	
	//	auto node = m_nodes[i];
	//	node->savePosition();
	//	//node->m_Jdot.clear();
	//}

	//// Integrate a small time step
	//yi = y0;
	//ydot0.segment(0, world->nr) = qdot;
	//yi += EPSILON * ydot0;
	//cout << "yi " << endl << yi << endl;


	//// Get new state
	//world->getJoint0()->scatterDofs(yi, world->nr);
	//world->getMuscle0()->update();

	//for (int t = 0; t < m_n_nodes; ++t) {
	//	m_nodes[t]->savePositionForJdot(); // the postion is saved in x_ss
	//}
	//// Compute Jdot for each node in the new time t0+EPS
	//for (int i = 0; i < m_n_bodies; i++) {
	//	// For each related body, perturb 
	//	yj = yi;
	//	int idx = m_bodies[i]->getJoint()->idxR;
	//	yj(idx) += EPSILON;
	//	
	//	world->getJoint0()->scatterDofs(yj, world->nr);
	//	world->getMuscle0()->update();

	//	for (int j = 0; j < m_n_nodes; ++j) {
	//		Vector3d diff = (m_nodes[j]->x - m_nodes[j]->x_ss) / EPSILON;
	//		m_nodes[j]->m_Jdot.col(idx) = diff ; // attention: the corresponding rigid body
	//	}
	//}

	//// Jdoti computed 
	//for (int j = 0; j < m_n_nodes; ++j) {
	//	auto node = m_nodes[j];

	//	//cout << "j:" << j << endl << "m_J" << endl << m_nodes[j]->m_J << endl << endl;
	//	//cout << "j:" << j << endl << "m_Jdot" << endl << m_nodes[j]->m_Jdot << endl << endl;
	//	MatrixXd diff = (node->m_J - node->m_Jdot) / EPSILON;
	//	//cout << "diff: " << endl <<  diff<< endl;
	//	m_nodes[j]->m_Jdot = diff;
	//}

	//// Restore the states
	//world->getJoint0()->scatterDofs(y0, world->nr);
	//world->getMuscle0()->update();


	//for (int i = 0; i < m_n_bodies; ++i) {
	//	// Compute Jdot for each node vector<3x2> 2 
	//	yi = y0;
	//	int idx = m_bodies[i]->getJoint()->idxR;
	//	yi(idx) += EPSILON;
	//	// before perturb
	//	world->getJoint0()->scatterDofs(yi, world->nr);
	//	world->getMuscle0()->update();
	//	for (int t = 0; t < m_n_nodes; ++t) {
	//		m_nodes[t]->savePositionForJdot();
	//	}
	//	// 

	//	for (int j = 0; j < m_n_bodies; j++) {
	//		// For each related body, perturb 
	//		yj = yi;
	//		int index = m_bodies[j]->getJoint()->idxR;
	//		yj(index) += EPSILON;
	//		//cout << "yj" << yj << endl;
	//		world->getJoint0()->scatterDofs(yj, world->nr);
	//		world->getMuscle0()->update();

	//		for (int k = 0; k < m_n_nodes; ++k) {
	//			Vector3d diff = (m_nodes[k]->x - m_nodes[k]->x_ss) / EPSILON;
	//			m_nodes[k]->m_Jdot[idx].col(index) = diff; // attention: the corresponding rigid body
	//		}
	//	}

	//	// Jdoti computed 
	//	// Restore the states
	//	world->getJoint0()->scatterDofs(y0, world->nr);
	//	world->getMuscle0()->update();
	//}

	//for (int t = 0; t < m_n_nodes; ++t) {
	//	for (int s = 0; s < m_n_bodies; ++s) { // each slice of Jdot
	//		m_nodes[t]->m_Jdot[s] = (m_nodes[t]->m_Jdot[s] - m_nodes[t]->m_J) / EPSILON;
	//		//cout << "node " << t << "slice " << s << endl << m_nodes[t]->m_Jdot[s] << endl << endl;
	//	}
	//}

	//MatrixXd Jdotqdot(3, m_n_bodies);
	//Jdotqdot.setZero();

	//for (int j = 0; j < m_n_nodes; ++j) {
	//	auto node = m_nodes[j];
	//	/*
	//	Jdotqdot.setZero();

	//	for (int k = 0; k < m_n_bodies; ++k) {
	//		int idx_k = m_bodies[k]->getJoint()->idxR;
	//		Jdotqdot += node->m_Jdot[idx_k] * qdot(idx_k);
	//	}*/
	//	// Jdotqdot for one node is done

	//	// compute JTMJDOTQDOT
	//	VectorXd fvec = node->m * node->m_J.transpose() * Matrix3d::Identity() * node->m_Jdot * qdot; // need to check the sign later
	//	
	//	f -= fvec;

	//}

}
