#pragma once
#include "JointUniversal.h"
#include "Body.h"
#include "SE3.h"
#include "Shape.h"

class JointUniversalXY : public JointUniversal {
public:
	JointUniversalXY() {};
	JointUniversalXY(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr):
	JointUniversal(body, parent)
	{

	}

	void update_() {
		double q0 = m_q(0);
		double q1 = m_q(1);
		double dq0 = m_qdot(0);
		double dq1 = m_qdot(1);

		double c0 = cos(q0);
		double c1 = cos(q1);
		double s0 = sin(q0);
		double s1 = sin(q1);

		m_Q = Matrix4d::Identity();
		Vector3d temp;
		temp << c1, s0*s1, -c0 * s1;
		m_Q.block<3, 1>(0, 0) = temp;
		temp << 0, c0, s0;
		m_Q.block<3, 1>(0, 1) = temp;
		temp << s1, -s0 * c1, c0 * c1;
		m_Q.block<3, 1>(0, 2) = temp;

		m_S(0, 0) = c1;
		m_S(2, 0) = s1;
		m_S(1, 1) = 1.0;
		m_Sdot(0, 0) = -s1 * dq1;
		m_Sdot(2, 0) = c1 * dq1;
	}

	virtual ~JointUniversalXY() {}

protected:


};