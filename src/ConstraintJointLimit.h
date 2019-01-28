#pragma once
// ConstraintJointLimit Inequality joint limit constraint
// Works only for revolute joints

#include "Constraint.h"

class Joint;

class ConstraintJointLimit : public Constraint 
{
public:
	ConstraintJointLimit();
	ConstraintJointLimit(std::shared_ptr<Joint> joint);
	void setLimits(double ql, double qu) { m_ql = ql; m_qu = qu; }

	std::shared_ptr<Joint> m_joint;	// Joint to apply the constriant on
	double m_ql;						// lower limit
	double m_qu;						// upper limit

protected:
	void computeJacIneqR_(Eigen::MatrixXd &Cr, Eigen::MatrixXd &Crdot, Eigen::VectorXd &cr, Eigen::VectorXd &crdot, Eigen::VectorXd &crddot);
	void ineqEventFcn_(std::vector<double> &value, std::vector<int> &isterminal, std::vector<int> &direction);
	void ineqProjPos_();
};