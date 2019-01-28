#pragma once
#include "Constraint.h"

class Joint;

class ConstraintPrescJoint : public Constraint, public std::enable_shared_from_this<ConstraintPrescJoint> 
{
public:
	ConstraintPrescJoint() {}
	ConstraintPrescJoint(std::shared_ptr<Joint> joint, Integrator vel);
	virtual ~ConstraintPrescJoint() {}
	void setActive() { activeER = true; }
	void setInactive() { activeER = false; }
	std::shared_ptr<Joint> m_joint;
	Eigen::VectorXd m_q;
	Eigen::VectorXd m_qdot;
	Eigen::VectorXd m_qddot;
	Integrator m_vel;
	std::shared_ptr<ConstraintPrescJoint> getSelf() { return shared_from_this(); }

protected:
	void computeJacEqR_(Eigen::MatrixXd &Gr, Eigen::MatrixXd &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void computeJacEqRSparse_(std::vector<T> &Gr, std::vector<T> &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void scatterForceEqR_();
	void init_();
private:

};