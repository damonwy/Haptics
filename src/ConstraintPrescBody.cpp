#include "rmpch.h"
#include "ConstraintPrescBody.h"
#include "Body.h"

using namespace std;
using namespace Eigen;

ConstraintPrescBody::ConstraintPrescBody(shared_ptr<Body> body, VectorXi prows, Integrator vel) :
	Constraint(prows.size(), 0, 0, 0), m_body(body), m_vel(vel), m_prows(prows)
{
	// set body ->presc = this;
	m_q.setZero();
	m_qdot.setZero();
	m_qddot.setZero();
	activeEM = false;

}
void ConstraintPrescBody::init_() {
	m_body->presc = getSelf();
}

void ConstraintPrescBody::computeJacEqM_(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
	int row = idxEM;
	int col = m_body->idxM;
	
	Matrix6d I6 = - Matrix6d::Identity();

	Gm.block(row, col, nconEM, 6) = I6(m_prows, Eigen::placeholders::all);

	if (m_vel == REDMAX_EULER) {
		gmdot.segment(row, nconEM) = m_qdot(m_prows);
	}
	else {
		gmdot.segment(row, nconEM) = m_body->phi(m_prows);
		gmddot.segment(row, nconEM) = m_qddot(m_prows);
	}
}

void ConstraintPrescBody::computeJacEqMSparse_(vector<T> &Gm, vector<T> &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
	int row = idxEM;
	int col = m_body->idxM;

	for (int i = 0; i < nconEM; ++i) {
		Gm.push_back(T(row + i, col + m_prows(i), -1.0));
	}

	if (m_vel == REDMAX_EULER) {
		gmdot.segment(row, nconEM) = m_qdot(m_prows);
	}
	else {
		gmdot.segment(row, nconEM) = m_body->phi(m_prows);
		gmddot.segment(row, nconEM) = m_qddot(m_prows);
	}
}

void ConstraintPrescBody::scatterForceEqM_() {
	m_body->wext_i = fcon;
}
