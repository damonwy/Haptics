#include "rmpch.h"
#include "Joint.h"
#include "Body.h"
#include "Constraint.h"
#include "ConstraintPrescJoint.h"

using namespace std;
using namespace Eigen;

Joint::Joint() {
	presc = nullptr;
}

Joint::Joint(shared_ptr<Body> body, int ndof, shared_ptr<Joint> parent) :
m_ndof(ndof),
m_body(body),
m_parent(parent)
{
	if (parent == nullptr) {
		m_name = "NULL-" + body->getName();
	}
	else {
		m_name = parent->getBody()->getName() + "-" + body->getName();
	}
	m_draw_radius = 1.0;
	m_q.resize(m_ndof);
	m_qdot.resize(m_ndof);
	m_qddot.resize(m_ndof);

	m_q.setZero();
	m_qdot.setZero();
	m_qddot.setZero();

	m_tau.resize(m_ndof);
	m_tau.setZero();

	m_tauCon.resize(m_ndof);
	m_tauCon.setZero();

	m_Kr = 0.0;
	m_Dr = 0.0;
	
	m_S.resize(6, ndof);
	m_S.setZero();
	m_Sdot.resize(6, ndof);
	m_Sdot.setZero();
	m_Q.setIdentity();
	
	m_I_j.setIdentity();
	V.setZero();
	Vdot.setZero();

	presc = nullptr;
}

void Joint::load(const string &RESOURCE_DIR, string joint_shape) {
	m_jointShape = make_shared<Shape>();
	m_jointShape->loadMesh(RESOURCE_DIR + joint_shape);

}

void Joint::init(int &nm, int &nr) {
	if (m_jointShape) {
		m_jointShape->init();
	}

	m_body->m_joint = getJoint();
	if (m_parent != nullptr) {
		m_parent->addChild(getJoint());
	}

	init_();

	countDofs(nm, nr);
}

void Joint::setJointTransform(Matrix4d E) {
	// Sets the transform of this joint wrt parent joint
	E_pj0 = E;
}

void Joint::update() {
	// Updates this joint and the attached body
	update_();
	// Transforms and adjoints
	E_pj.noalias() = E_pj0 * m_Q;

	E_jp = SE3::inverse(E_pj);
	Ad_jp = SE3::adjoint(E_jp);

	Matrix4d E_wp;

	if (m_parent == nullptr) {
		E_wp.setIdentity();
	}
	else {
		E_wp = m_parent->E_wj;
	}

	E_wj.noalias() = E_wp * E_pj;
	// Joint velocity
	V.noalias() = m_S * m_qdot;

	if (m_parent != nullptr) {
		// Add parent velocity
		V.noalias() += Ad_jp * m_parent->V;
	}

	// Update attached body
	m_body->update();
	if (next != nullptr) {
		next->update();
	}

}

void Joint::countDofs(int &nm, int &nr) {
	// Counts reduced DOFs
	idxR = countR(nr, m_ndof);
	m_body->countDofs(nm);
	m_q0 = m_q;
}

void Joint::countHRDofs(int &nR) {
	idxHR = countR(nR, m_ndof);
}



int Joint::countR(int &nr, int data) {
	nr = nr + data;
	return (nr - data);
}

void Joint::computeJacobian(MatrixXd &J, MatrixXd &Jdot) {
	// Computes the redmax Jacobian
	Matrix6d Ad_ij = m_body->Ad_ij;
	J.block(m_body->idxM, idxR, 6, m_ndof).noalias() = Ad_ij * m_S;
	Jdot.block(m_body->idxM, idxR, 6, m_ndof).noalias() = Ad_ij * m_Sdot;

	// Loop through all ancestors
	auto jointA = m_parent;
	while (jointA != nullptr) {
		int idxM_P = m_parent->getBody()->idxM;
		Matrix6d Ad_ip = m_body->Ad_ip;
		Matrix6d Ad_iw = m_body->Ad_iw;
		Matrix6d Ad_wp = m_parent->getBody()->Ad_wi;
		Matrix6d Addot_wi = m_body->Addot_wi;
		Matrix6d Addot_wp = m_parent->getBody()->Addot_wi;
		Matrix6d Addot_ip = -Ad_iw * (Addot_wi * Ad_iw * Ad_wp - Addot_wp);
		J.block(m_body->idxM, jointA->idxR, 6, jointA->m_ndof) = Ad_ip * J.block(idxM_P, jointA->idxR, 6, jointA->m_ndof);
		Jdot.block(m_body->idxM, jointA->idxR, 6, jointA->m_ndof) = Ad_ip * Jdot.block(idxM_P, jointA->idxR, 6, jointA->m_ndof) + Addot_ip * J.block(idxM_P, jointA->idxR, 6, jointA->m_ndof);
		jointA = jointA->getParent();
	}


	if (next != nullptr) {
		next->computeJacobian(J, Jdot);
	}
}

void Joint::computeHyperReducedJacobian(MatrixXd &JrR, MatrixXd &JrR_select) {
	// Computes the chain Hyper Reduced Jacobian JmR
	computeHyperReducedJacobian_(JrR, JrR_select);
	if (next != nullptr) {
		next->computeHyperReducedJacobian(JrR, JrR_select);
	}
}

void Joint::computeHyperReducedJacobian_(MatrixXd &JrR, MatrixXd &JrR_select) {
	for (int i = 0; i < m_ndof; ++i) {
		JrR(idxR + i, idxHR + i) = 1.0;
		JrR_select(idxR + i, idxHR + i) = 1.0;
	}
}


void Joint::computeInertia() {
	double m = m_body->I_i(3);

	Matrix3d R, pBrac, Ic;
	Vector3d p;
	SE3::EToRp(m_body->E_ji, R, p);

	pBrac = SE3::bracket3(p);
	Ic = (m_body->I_i.segment<3>(0)).asDiagonal();

	m_I_j.block<3, 3>(0, 0).noalias() = R * Ic *R.transpose() + m *(pBrac.transpose()*pBrac);
	m_I_j.block<3, 3>(0, 3) = m * pBrac;
	m_I_j.block<3, 3>(3, 0).noalias() = m * pBrac.transpose();
	m_I_j.block<3, 3>(3, 3) = m * Matrix3d::Identity();
}

void Joint::reparam() {
	reparam_();
	if (next != nullptr) {
		next->reparam();
	}
}

void Joint::computeForceStiffness(VectorXd &fr, MatrixXd &Kr) {
	// Computes joint stiffness force vector and matrix
	if (presc == nullptr) {
		int row = this->idxR;
		// Add the joint torque here rather than having a separate function
		fr.segment(row, m_ndof).noalias() += m_tau - m_Kr * (m_q - m_q0);
		MatrixXd I(m_ndof, m_ndof);
		I.setIdentity();
		Kr.block(row, row, m_ndof, m_ndof) -= m_Kr * I;
	}

	if (next != nullptr) {
		next->computeForceStiffness(fr, Kr);
	}
}

void Joint::computeForceStiffnessSparse(VectorXd &fr, vector<T> &Kr_) {
	// Computes joint stiffness force vector and matrix
	if (presc == nullptr) {
		int row = this->idxR;
		// Add the joint torque here rather than having a separate function
		fr.segment(row, m_ndof) += m_tau - m_Kr * m_q;

		for (int i = 0; i < m_ndof; ++i) {
			Kr_.push_back(T(row + i, row + i, - m_Kr));
		}
	}

	if (next != nullptr) {
		next->computeForceStiffnessSparse(fr, Kr_);
	}
}

void Joint::computeForceDamping(VectorXd &fr, MatrixXd &Dr) {
	// Computes joint damping force vector and matrix
	if (presc == nullptr) {
		int row = this->idxR;
		fr.segment(row, m_ndof).noalias() -= m_Dr * m_qdot;
		MatrixXd I(m_ndof, m_ndof);
		I.setIdentity();
		Dr.block(row, row, m_ndof, m_ndof) += m_Dr * I;
	}

	if (next != nullptr) {
		next->computeForceDamping(fr, Dr);
	}
}

void Joint::computeForceDampingSparse(VectorXd &fr, vector<T> &Dr_) {
	if (presc == nullptr) {
		int row = this->idxR;
		fr.segment(row, m_ndof) -= m_Dr * m_qdot;

		for (int i = 0; i < m_ndof; ++i) {
			Dr_.push_back(T(row + i, row + i, m_Dr));
		}
	}

	if (next != nullptr) {
		next->computeForceDampingSparse(fr, Dr_);
	}
}

VectorXd Joint::computerJacTransProd(VectorXd y, VectorXd x, int nr) {
	// Computes x = J'*y
	// x (nr, 1)
	VectorXd yi = y.segment(m_body->idxM, 6);
	for (int k = 0; k < (int)m_children.size(); k++) {
		yi = yi + m_children[k]->getAlpha();
	}
	m_alpha = m_body->Ad_ip.transpose() * yi;
	x.segment(idxR, m_ndof) = (m_body->Ad_ij * m_S).transpose() * yi;

	if (prev != nullptr) {
		x = prev->computerJacTransProd(y, x, nr);
	}
	return x;
}

void Joint::computeEnergies(Vector3d grav, Energy &ener) {
	// Computes kinetic and potential energies
	m_body->computeEnergies(grav, ener);
	VectorXd dq = m_q - m_q0;
	ener.V += 0.5 * m_Kr * dq.dot(dq);

	if (next != nullptr) {
		next->computeEnergies(grav, ener);
	}
}

void Joint::gatherDofs(VectorXd &y, int nr) {
	// Gathers q and qdot into y
	y.segment(idxR, m_ndof) = m_q;
	y.segment(nr + idxR, m_ndof) = m_qdot;
	if (next != nullptr) {
		next->gatherDofs(y, nr);
	}
}

void Joint::gatherDDofs(VectorXd &ydot, int nr) {
	// Gathers qdot and qddot into ydot
	ydot.segment(idxR, m_ndof) = m_qdot;
	ydot.segment(nr + idxR, m_ndof) = m_qddot;
	if (next != nullptr) {
		next->gatherDDofs(ydot, nr);
	}
}

void Joint::scatterDofs(VectorXd y, int nr) {
	// Scatters q and qdot from y
	scatterDofsNoUpdate(y, nr);
	update();
}

void Joint::scatterDDofs(VectorXd ydot, int nr) {
	// Scatters qdot and qddot from ydot
	m_qdot.segment(0, m_ndof) = ydot.segment(idxR, m_ndof);
	m_qddot.segment(0, m_ndof) = ydot.segment(nr + idxR, m_ndof);
	if (next != nullptr) {
		next->scatterDDofs(ydot, nr);
	}
}

void Joint::scatterDofsNoUpdate(VectorXd y, int nr) {
	// Helper function to scatter without updating
	m_q.segment(0, m_ndof) = y.segment(idxR, m_ndof);
	m_qdot.segment(0, m_ndof) = y.segment(nr + idxR, m_ndof);
	if (next != nullptr) {
		next->scatterDofsNoUpdate(y, nr);
	}
}

void Joint::scatterTauCon(VectorXd tauc) {
	// Scatters constraint force
	m_tauCon = tauc.segment(idxR, m_ndof);
	if (next != nullptr) {
		next->scatterTauCon(tauc);
	}
}

void Joint::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {

	progSimple->bind();
	MV->pushMatrix();
	MV->multMatrix(eigen_to_glm(E_wj));
	glUniformMatrix4fv(progSimple->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	
	glLineWidth(5);
	glBegin(GL_LINES);
	// X axis
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(3.0f, 0.0f, 0.0f);

	// Y axis
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 3.0f, 0.0f);

	// Z axis
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 3.0f);

	glEnd();
	MV->popMatrix();
	progSimple->unbind();

	draw_(MV, prog, progSimple, P);

	if (next != nullptr) {
		next->draw(MV, prog, progSimple, P);
	}

}

void Joint::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
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
		if (presc->activeER) {
			alpha = 2.0;
		}
		MV->scale(static_cast<float>(alpha * m_draw_radius));
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		m_jointShape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();
}
