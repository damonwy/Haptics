#pragma once

#ifndef REDUCEDCOORD_SRC_JOINTSPHERICALEXP_H_
#define REDUCEDCOORD_SRC_JOINTSPHERICALEXP_H_

#include "Joint.h"
#include "SE3.h"
#include "Body.h"
#include "MatrixStack.h"
#include "Program.h"

class JointSphericalExp : public Joint {

public:
	JointSphericalExp() {}
	JointSphericalExp(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr):
	Joint(body, 3, parent)
	{
		m_radius = 1.0;
	}

	void load(const std::string &RESOURCE_DIR, std::string joint_shape) {

	}

	virtual ~JointSphericalExp() {}
	void update_() {
		Vector3d r = m_q.topRows(3);
		Vector3d rdot = m_qdot.topRows(3);
		m_Q = Matrix4d::Identity();
		Matrix3d R = SE3::exp(r);
		m_Q.block<3, 3>(0, 0) = R;
		m_S.setZero();
		m_Sdot.setZero();
		double rr = r.dot(r);
		if (rr < 1.0e-9) {
			m_S.block<3, 3>(0, 0).setIdentity();
		}
		else {
			Matrix3d Rdot = Matrix3d::Zero();
			Matrix3d rbrac = SE3::bracket3(r);
			Matrix3d IR = Matrix3d::Identity() - R;
			double d = 1.0 / rr;
			for (int i = 0; i < 3; ++i) {
				Matrix3d Bi = r(i)*rbrac;
				Matrix3d Ci = SE3::bracket3(rbrac*IR.col(i));
				Matrix3d Ai = (Bi + Ci)*d;
				Matrix3d dRdri = Ai * R;
				m_S.block<3, 1>(0, i) = SE3::unbracket3(R.transpose()*dRdri);
				Rdot += dRdri *rdot(i);
			}

			Matrix3d rdotbrac = SE3::bracket3(rdot);
			double ddot = -2.0 / (rr * rr) * (r.dot(rdot));

			for (int i = 0; i < 3; ++i) {
				Matrix3d Bi = r(i)*rbrac;
				Matrix3d Ci = SE3::bracket3(rbrac*IR.col(i));
				Matrix3d Ai = (Bi + Ci)*d;
				Matrix3d Bidot = rdot(i) * rbrac + r(i) * rdotbrac;
				Matrix3d Cidot = SE3::bracket3(rdotbrac * IR.col(i) - rbrac * Rdot.col(i));
				Matrix3d Aidot = (Bidot + Cidot)*d + (Bi + Ci)*ddot;
				Matrix3d RdotAiR = Rdot.transpose() * Ai * R;

				m_Sdot.block<3, 1>(0, i) = SE3::unbracket3(RdotAiR + R.transpose()* Aidot * R - RdotAiR.transpose());
			}
		}

	}

	void reparam_() {
		Vector3d q0 = m_q.topRows(3);
		Eigen::VectorXd q1 = q0;
		bool flag = SE3::reparam(q1);

		if (flag) {
			Vector3d qdot0 = m_qdot.topRows(3);
			Matrix3d S0 = m_S.block<3, 3>(0, 0);
			// Update to compute new S (ignore new Sdot)
			m_q.topRows(3) = q1;
			update_();
			Matrix3d S1 = m_S.block<3, 3>(0, 0);

			// Update again to compute new Sdot
			m_qdot.topRows(3) = S1.partialPivLu().solve(S0*qdot0);
			update_();
		}
	}

	void setGeometry(double radius) {
		m_radius = radius;
	}

protected:

	void draw_(std::shared_ptr<MatrixStack> MV,
		const std::shared_ptr<Program> prog,
		const std::shared_ptr<Program> prog2,
		std::shared_ptr<MatrixStack> P) const {
	}

	double m_radius;	// Radius for contact

};



#endif // REDUCEDCOORD_SRC_JOINTSPHERICALEXP_H_