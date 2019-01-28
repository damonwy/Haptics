#pragma once
#include "Constraint.h"

class Body;

class ConstraintPrescBodyAttachPoint : public Constraint, public std::enable_shared_from_this<ConstraintPrescBodyAttachPoint> {
public:
	ConstraintPrescBodyAttachPoint() {}
	ConstraintPrescBodyAttachPoint(std::shared_ptr<Body> body, Vector3d r, Eigen::VectorXi prows, Integrator vel);
	std::shared_ptr<ConstraintPrescBodyAttachPoint> getSelf() { return shared_from_this(); }
	virtual ~ConstraintPrescBodyAttachPoint() {}
	Vector3d m_q;
	Vector3d m_qdot;
	Vector3d m_qddot;
	void setActive() { activeEM = true; }
	void setInactive() { activeEM = false; }

protected:
	void computeJacEqM_(Eigen::MatrixXd &Gr, Eigen::MatrixXd &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void computeJacEqMSparse_(std::vector<T> &Gr, std::vector<T> &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void scatterForceEqM_();
	void init_();

private:
	std::shared_ptr<Body> m_body;
	Vector3d m_r; // local position of the point
	Integrator m_vel;
	Eigen::VectorXi m_prows;
	Matrix3x6d m_gamma;
};