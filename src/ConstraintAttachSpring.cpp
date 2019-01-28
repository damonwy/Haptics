#include "rmpch.h"
#include "ConstraintAttachSpring.h"

#include "Deformable.h"
#include "DeformableSpring.h"
#include "Node.h"
#include "Body.h"

using namespace std;
using namespace Eigen;

ConstraintAttachSpring::ConstraintAttachSpring() {

}

ConstraintAttachSpring::ConstraintAttachSpring(shared_ptr<Deformable> spring) :
Constraint(6, 0, 0, 0), m_spring(spring)
{

}

void ConstraintAttachSpring::computeJacEqM_(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
	int row0 = idxEM;
	int row1 = idxEM + 3;
	int col0S = m_spring->m_nodes[0]->idxM;
	int col1S = m_spring->m_nodes[m_spring->m_nodes.size()-1]->idxM;
	auto body0 = m_spring->m_body0;
	auto body1 = m_spring->m_body1;

	Matrix4d E0, E1;
	int col0B, col1B;
	if (body0 == nullptr) {
		E0 = Matrix4d::Identity();
		
	}
	else {
		E0 = body0->E_wi;
		col0B = body0->idxM;
	}

	if (body1 == nullptr) {
		E1 = Matrix4d::Identity();

	}
	else {
		E1 = body1->E_wi;
		col1B = body1->idxM;
	}

	Matrix3d R0 = E0.block<3, 3>(0, 0);
	Matrix3d R1 = E1.block<3, 3>(0, 0);

	Matrix3x6d G0 = SE3::gamma(m_spring->m_r0);
	Matrix3x6d G1 = SE3::gamma(m_spring->m_r1);

	Matrix3d W0, W1;
	if (body0 != nullptr) {
		W0 = SE3::bracket3(body0->phi.segment<3>(0));
		Gm.block<3, 6>(row0, col0B) = R0 * G0;
		Gmdot.block<3, 6>(row0, col0B) = R0 * W0 * G0;
	}

	if (body1 != nullptr) {
		W1 = SE3::bracket3(body1->phi.segment<3>(0));
		Gm.block<3, 6>(row1, col1B) = R1 * G1;
		Gmdot.block<3, 6>(row1, col1B) = R1 * W1 * G1;
	}

	Gm.block<3, 3>(row0, col0S) = -Matrix3d::Identity();
	Gm.block<3, 3>(row1, col1S) = -Matrix3d::Identity();

	Vector4d tem00;
	tem00.segment<3>(0) = m_spring->m_r0;
	tem00(3) = 1.0;
	Vector4d tem01;
	tem01.segment<3>(0) = m_spring->m_nodes[0]->x;
	tem01(3) = 1.0;


	Vector4d tem10;
	tem10.segment<3>(0) = m_spring->m_r1;
	tem10(3) = 1.0;
	Vector4d tem11;
	tem11.segment<3>(0) = m_spring->m_nodes[m_spring->m_nodes.size()-1]->x;
	tem11(3) = 1.0;

	Vector4d gm0 = E0 * tem00 - tem01;
	Vector4d gm1 = E1 * tem10 - tem11;

	gm.segment<3>(row0) = gm0.segment<3>(0);
	gm.segment<3>(row1) = gm1.segment<3>(0);
}

void ConstraintAttachSpring::computeJacEqMSparse_(vector<T> &Gm, vector<T> &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
	int row0 = idxEM;
	int row1 = idxEM + 3;
	int col0S = m_spring->m_nodes[0]->idxM;
	int col1S = m_spring->m_nodes[m_spring->m_nodes.size() - 1]->idxM;
	auto body0 = m_spring->m_body0;
	auto body1 = m_spring->m_body1;

	Matrix4d E0, E1;
	int col0B, col1B;
	if (body0 == nullptr) {
		E0 = Matrix4d::Identity();

	}
	else {
		E0 = body0->E_wi;
		col0B = body0->idxM;
	}

	if (body1 == nullptr) {
		E1 = Matrix4d::Identity();

	}
	else {
		E1 = body1->E_wi;
		col1B = body1->idxM;
	}

	Matrix3d R0 = E0.block<3, 3>(0, 0);
	Matrix3d R1 = E1.block<3, 3>(0, 0);

	Matrix3x6d G0 = SE3::gamma(m_spring->m_r0);
	Matrix3x6d G1 = SE3::gamma(m_spring->m_r1);

	Matrix3d W0, W1;
	if (body0 != nullptr) {
		W0 = SE3::bracket3(body0->phi.segment<3>(0));

		Matrix3x6d Gm_temp = R0 * G0;
		Matrix3x6d Gmdot_temp = R0 * W0 * G0;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 6; j++) {
				Gm.push_back(T(row0 + i, col0B + j, Gm_temp(i, j)));
				Gmdot.push_back(T(row0 + i, col0B + j, Gmdot_temp(i, j)));
			}
		}
	}

	if (body1 != nullptr) {
		W1 = SE3::bracket3(body1->phi.segment<3>(0));
		Matrix3x6d Gm_temp = R1 * G1;
		Matrix3x6d Gmdot_temp = R1 * W1 * G1;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 6; j++) {
				Gm.push_back(T(row1 + i, col1B + j, Gm_temp(i, j)));
				Gmdot.push_back(T(row1 + i, col1B + j, Gmdot_temp(i, j)));
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		Gm.push_back(T(row0 + i, col0S + i, -1.0));
		Gm.push_back(T(row1 + i, col1S + i, -1.0));
	}

	Vector4d tem00;
	tem00.segment<3>(0) = m_spring->m_r0;
	tem00(3) = 1.0;
	Vector4d tem01;
	tem01.segment<3>(0) = m_spring->m_nodes[0]->x;
	tem01(3) = 1.0;

	Vector4d tem10;
	tem10.segment<3>(0) = m_spring->m_r1;
	tem10(3) = 1.0;
	Vector4d tem11;
	tem11.segment<3>(0) = m_spring->m_nodes[m_spring->m_nodes.size() - 1]->x;
	tem11(3) = 1.0;

	Vector4d gm0 = E0 * tem00 - tem01;
	Vector4d gm1 = E1 * tem10 - tem11;

	gm.segment<3>(row0) = gm0.segment<3>(0);
	gm.segment<3>(row1) = gm1.segment<3>(0);
}
