#include "rmpch.h"
#include "ConstraintPrescBodyAttachPoint.h"
#include "Body.h"

using namespace std;
using namespace Eigen;

ConstraintPrescBodyAttachPoint::ConstraintPrescBodyAttachPoint(shared_ptr<Body> body, Vector3d r, VectorXi prows, Integrator vel) :
	Constraint(prows.size(), 0, 0, 0), m_body(body), m_r(r), m_prows(prows), m_vel(vel)
{
	// set body ->presc = this;
	m_q.setZero();
	m_qdot.setZero();
	m_qddot.setZero();
	activeEM = false;
	m_gamma.setZero();
	m_gamma.leftCols(3) = SE3::bracket3(m_r).transpose();
	m_gamma.rightCols(3) = Matrix3d::Identity();

}
void ConstraintPrescBodyAttachPoint::init_() {
	m_body->m_presc_attach_points.push_back(getSelf());
}

void ConstraintPrescBodyAttachPoint::computeJacEqM_(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
	int row = idxEM;
	int col = m_body->idxM;
	Matrix3d R = m_body->E_wi.topLeftCorner(3, 3);
	Matrix3x6d Cons = -R * m_gamma;

	Gm.block(row, col, nconEM, 6) = Cons(m_prows, Eigen::placeholders::all);

	if (m_vel == REDMAX_EULER) {
		gmdot.segment(row, nconEM) = m_qdot(m_prows);
	}
	else {
		gmdot.segment(row, nconEM) = m_body->phi(m_prows);
		gmddot.segment(row, nconEM) = m_qddot(m_prows);
	}
}

void ConstraintPrescBodyAttachPoint::computeJacEqMSparse_(vector<T> &Gm, vector<T> &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
	int row = idxEM;
	int col = m_body->idxM;
	Matrix3d R = m_body->E_wi.topLeftCorner(3, 3);
	Matrix3x6d Cons = -R * m_gamma;
	//cout << "vel now : " << Cons * m_body->phi << endl;
	//Vector3d qdot1 = -Cons * m_body->phi; // drift?
	//Vector3d diff = m_qdot - qdot1;
	//cout << "diff " << diff << endl;

	for (int i = 0; i < nconEM; ++i) {
		Vector6d con = Cons.row(m_prows(i));

		for (int j = 0; j < 6; ++j) {
			Gm.push_back(T(row + i, col + j, con(j)));
		}
	}

	if (m_vel == REDMAX_EULER) {
		gmdot.segment(row, nconEM) = m_qdot(m_prows) ;
	}
	else {
		gmdot.segment(row, nconEM) = m_body->phi(m_prows);
		gmddot.segment(row, nconEM) = m_qddot(m_prows);
	}
}

void ConstraintPrescBodyAttachPoint::scatterForceEqM_() {
	//m_body->wext_i = fcon;
}