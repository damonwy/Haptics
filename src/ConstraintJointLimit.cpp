#include "rmpch.h"
#include "ConstraintJointLimit.h"
#include "Joint.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

ConstraintJointLimit::ConstraintJointLimit() {
}

ConstraintJointLimit::ConstraintJointLimit(std::shared_ptr<Joint> joint):
Constraint(0, 0, 0, 1),
m_joint(joint)
{
	m_name =  m_joint->getName() + "-LIMIT";

}

void ConstraintJointLimit::computeJacIneqR_(Eigen::MatrixXd &Cr, Eigen::MatrixXd &Crdot, Eigen::VectorXd &cr, Eigen::VectorXd &crdot, Eigen::VectorXd &crddot) {
	int row = idxIR;
	int col = m_joint->idxR;
	//idxQ = col;
	nQ = m_joint->m_ndof;
	idxQ.resize(nQ, 1);
	int temp = col;
	for (int i = 0; i < nQ; i++) {	
		idxQ(i) = temp;
		temp++;
	}

	if (m_joint->m_q(0) <= m_ql) {
		Cr.block(row, col, nconIR, m_joint->m_ndof).setOnes();
		Cr.block(row, col, nconIR, m_joint->m_ndof) *= -1;
		cr(row) = m_ql - m_joint->m_q(0);
		activeR = true;
	}
	else if (m_joint->m_q(0) >= m_qu) {
		Cr.block(row, col, nconIR, m_joint->m_ndof).setOnes();
		cr(row) = m_qu - m_joint->m_q(0);
		activeR = true;
	}
	else {
		activeR = false;
	}
}

void ConstraintJointLimit::ineqEventFcn_(vector<double> &value, vector<int> &isterminal, vector<int> &direction) {
	double q = m_joint->m_q(0);
	value.push_back(q - m_ql);
	isterminal.push_back(1);
	direction.push_back(-1);
	value.push_back(q - m_qu);
	isterminal.push_back(1);
	direction.push_back(1);
}

void ConstraintJointLimit::ineqProjPos_() {
	if (m_joint->m_q(0) <= m_ql) {
		m_joint->m_q(0) = m_ql;
	}
	else if (m_joint->m_q(0) >= m_qu) {
		m_joint->m_q(0) = m_qu;
	}

}