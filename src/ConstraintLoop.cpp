#include "rmpch.h"
#include "ConstraintLoop.h"
#include "Joint.h"
#include "Body.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

ConstraintLoop::ConstraintLoop() {
}

ConstraintLoop::ConstraintLoop(shared_ptr<Body> bodyA, shared_ptr<Body> bodyB) :
	Constraint(2, 0, 0, 0),
	m_bodyA(bodyA), m_bodyB(bodyB)
{
	m_name = m_bodyA->getName() + "-" + m_bodyB->getName();

}

void ConstraintLoop::computeJacEqM_(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
	int row = idxEM;
	Matrix4d E_wa = m_bodyA->E_wi;
	Matrix4d E_wb = m_bodyB->E_wi;
	Matrix3d R_wa = E_wa.block<3, 3>(0, 0);
	Matrix3d R_wb = E_wb.block<3, 3>(0, 0);

	// Get two directions orthonormal to A's hinge axis
	auto jointA = m_bodyA->getJoint();
	Vector3d v0 = R_wa * jointA->m_axis;

	Vector3d v0_abs = v0.cwiseAbs();
	// Get Location of minimum
	VectorXd::Index min_i;
	Vector3d v1; 
	v1.setZero();
	v1(min_i) = 1.0;
	
	Vector3d v2 = v0.cross(v1);
	v2.normalized();

	v1 = v2.cross(v0);
	v1.normalized();
	MatrixXd v12(3, 2);
	v12 << v1, v2;
	
	Matrix3x6d GammaA = SE3::gamma(m_xA);
	Matrix3x6d GammaB = SE3::gamma(m_xB);

	Matrix3d waBrac = SE3::bracket3(m_bodyA->phi.segment<3>(0));
	Matrix3d wbBrac = SE3::bracket3(m_bodyB->phi.segment<3>(0));
	// idxQ - [colsA, colsB];
	
	int colA = m_bodyA->idxM;
	int colB = m_bodyB->idxM;
	idxQ.resize(6, 2);
	idxQ.col(0) << colA, colA + 1, colA + 2, colA + 3, colA + 4, colA + 5;
	idxQ.col(1) << colB, colB + 1, colB + 2, colB + 3, colB + 4, colB + 5;

	Gm.block(row, colA, nconEM, 6) = v12.transpose() * R_wa * GammaA;
	Gm.block(row, colB, nconEM, 6) = -v12.transpose() * R_wb * GammaB;

	Gmdot.block(row, colA, nconEM, 6) = v12.transpose() * R_wa * waBrac * GammaA;
	Gmdot.block(row, colB, nconEM, 6) = -v12.transpose() * R_wb * wbBrac * GammaB;

	Vector4d temp0, temp1;
	temp0 << m_xA, 1.0;
	temp1 << m_xB, 1.0;

	Vector4d dx = E_wa * temp0 - E_wb * temp1;
	gm.segment<2>(row) = v12.transpose() * dx.segment<3>(0);
}
