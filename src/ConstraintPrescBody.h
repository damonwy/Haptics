#pragma once
#include "Constraint.h"

class Body;

class ConstraintPrescBody : public Constraint, public std::enable_shared_from_this<ConstraintPrescBody> {
public:
	ConstraintPrescBody() {}
	ConstraintPrescBody(std::shared_ptr<Body> body, Eigen::VectorXi prows, Integrator vel);	
	std::shared_ptr<ConstraintPrescBody> getSelf() { return shared_from_this(); }
	virtual ~ConstraintPrescBody() {}
	Vector6d m_q;
	Vector6d m_qdot;
	Vector6d m_qddot;	
	void setActive() { activeEM = true; }
	void setInactive() { activeEM = false; }

protected:
	void computeJacEqM_(Eigen::MatrixXd &Gr, Eigen::MatrixXd &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void computeJacEqMSparse_(std::vector<T> &Gr, std::vector<T> &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void scatterForceEqM_();
	void init_();

private:
	std::shared_ptr<Body> m_body;

	Integrator m_vel;
	Eigen::VectorXi m_prows;
};