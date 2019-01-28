#pragma once
// ConstraintAttachSpring

#include "Constraint.h"

class Deformable;
class DeformableSpring;

class ConstraintAttachSpring : public Constraint
{
public:
	ConstraintAttachSpring();
	ConstraintAttachSpring(std::shared_ptr<Deformable> spring);
	std::shared_ptr<Deformable> m_spring;

protected:
	void computeJacEqM_(Eigen::MatrixXd &Gm, Eigen::MatrixXd &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);
	void computeJacEqMSparse_(std::vector<T> &Gm, std::vector<T> &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);

};