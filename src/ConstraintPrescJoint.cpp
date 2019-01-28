#include "rmpch.h"
#include "ConstraintPrescJoint.h"
#include "Joint.h"

using namespace std;
using namespace Eigen;

ConstraintPrescJoint::ConstraintPrescJoint(shared_ptr<Joint> joint, Integrator vel):
Constraint(0, joint->m_ndof, 0, 0), m_joint(joint), m_vel(vel)
{
	m_q.resize(joint->m_ndof);
	m_qdot.resize(joint->m_ndof);
	m_qddot.resize(joint->m_ndof);
	m_q.setZero();
	m_qdot.setZero();
	m_qddot.setZero();
	activeER = false;
}

void ConstraintPrescJoint::init_() {
	m_joint->presc = getSelf();
}

void ConstraintPrescJoint::computeJacEqR_(MatrixXd &Gr, MatrixXd &Grdot, VectorXd &gr, VectorXd &grdot, VectorXd &grddot) {
	int row = idxER;
	int col = m_joint->idxR;
	nidxQ = m_joint->m_ndof;
	MatrixXd In(m_joint->m_ndof, m_joint->m_ndof);
	In.setIdentity();
	In *= -1.0;

	Gr.block(row, col, nconER, nidxQ) = In;
	gr.segment(row, nconER) = m_q - m_joint->m_q;
	if (m_vel == REDMAX_EULER) {
		grdot.segment(row, nconER) = m_qdot;
	}
	else {
		grdot.segment(row, nconER) = m_joint->m_qdot;
		grddot.segment(row, nconER) = m_qddot;
	}
}

void ConstraintPrescJoint::computeJacEqRSparse_(vector<T> &Gr, vector<T> &Grdot, VectorXd &gr, VectorXd &grdot, VectorXd &grddot) {
	int row = idxER;
	int col = m_joint->idxR;
	nidxQ = m_joint->m_ndof;
	MatrixXd In(m_joint->m_ndof, m_joint->m_ndof);
	In.setIdentity();
	In *= -1.0;

	//Gr.block(row, col, nconER, nidxQ) = In;
	for (int i = 0; i < m_joint->m_ndof; ++i) {
		Gr.push_back(T(row + i, col + i, -1.0));
	}

	gr.segment(row, nconER) = m_q - m_joint->m_q;
	
	if (m_vel == REDMAX_EULER) {
		grdot.segment(row, nconER) = m_qdot;
	}
	else {
		grdot.segment(row, nconER) = m_joint->m_qdot;
		grddot.segment(row, nconER) = m_qddot;
	}
}


void ConstraintPrescJoint::scatterForceEqR_() {
	m_joint->m_tau = fcon;
}